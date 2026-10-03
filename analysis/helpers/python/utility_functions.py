# -*- coding: utf-8 -*-
"""
Created on Sat Jun 28 13:22:31 2025

@author: Fumiya
"""


import numpy as np
import matplotlib.pyplot as plt  # for plotting
import os
import time
from matplotlib.patches import Wedge


def BlockCat(data: np.ndarray, dim2squish: int, dim2concat: int) -> np.ndarray:
	"""
	Split `data` along axis `dim2squish` into single‐slice blocks,
	then concatenate those blocks along axis `dim2concat`.

	Parameters
	----------
	data : np.ndarray
		The array to be rearranged.
	dim2squish : int
		The axis along which to split into blocks of size 1.
	dim2concat : int
		The axis along which to concatenate all those blocks.

	Returns
	-------
	result : np.ndarray
		The resulting array.
	"""
	# normalize negative axes
	ndim = data.ndim
	i = dim2squish if dim2squish >= 0 else ndim + dim2squish
	j = dim2concat if dim2concat >= 0 else ndim + dim2concat
	if not (0 <= i < ndim and 0 <= j < ndim):
		raise np.AxisError(f"axes out of range (data.ndim={ndim}): "
						   f"dim2squish={dim2squish}, dim2concat={dim2concat}")

	# number of slices along the split axis
	n_slices = data.shape[i]

	# split into n_slices blocks (each of size 1 along axis i)
	blocks = np.split(data, n_slices, axis=i)

	# concatenate all the blocks along axis j
	return np.concatenate(blocks, axis=j)

def invBlockCat(
		concatenated: np.ndarray,
		original:     np.ndarray,
		dim2squish:   int,
		dim2concat:   int
	) -> np.ndarray:
	"""
	Reverse of BlockCat: given the array after BlockCat (concatenated),
	and the original array (to recover its shape), split along dim2concat
	and re‐assemble along dim2squish.

	Parameters
	----------
	concatenated : np.ndarray
		The result of BlockCat(original, dim2squish, dim2concat).
	original : np.ndarray
		The array that was fed into BlockCat.
	dim2squish : int
		The axis along which original was split into unit‐sized blocks.
	dim2concat : int
		The axis along which those blocks were concatenated.

	Returns
	-------
	np.ndarray
		The reconstructed original array.
	"""
	# normalize negative axes
	ndim = original.ndim
	i = dim2squish if dim2squish >= 0 else ndim + dim2squish
	j = dim2concat if dim2concat >= 0 else ndim + dim2concat
	if not (0 <= i < ndim and 0 <= j < ndim):
		raise np.AxisError(f"axes out of range (ndim={ndim}): "
						   f"dim2squish={dim2squish}, dim2concat={dim2concat}")

	# how many unit‐blocks there were
	n_slices = original.shape[i]

	# sanity check: concatenated.shape[j] should equal
	# original.shape[j] * n_slices  (except when i==j)
	expected = (original.shape[j] * n_slices
				if i != j else
				original.shape[j])
	if concatenated.shape[j] != expected:
		raise ValueError(
			f"Unexpected concat‐axis size: got {concatenated.shape[j]}, "
			f"expected {expected}"
		)

	# split back into the n_slices blocks along axis j...
	blocks = np.split(concatenated, n_slices, axis=j)
	# ...and glue them back along axis i
	return np.concatenate(blocks, axis=i)

def multiBlockCat(data: np.ndarray,
				  dim2squish: int | list[int],
				  dim2concat: int) -> np.ndarray:
	"""
	Repeatedly apply BlockCat to `data`.

	Parameters
	----------
	data : np.ndarray
		The array to reshape.
	dim2squish : int or list of int
		Axis or axes along which to split into single‐slice blocks.
		If a list, BlockCat is applied in the given order.
	dim2concat : int
		The axis along which to concatenate after each split.

	Returns
	-------
	result : np.ndarray
		The reshaped array after all splits and concatenations.
	"""
	# ensure we have a list of axes
	axes = [dim2squish] if isinstance(dim2squish, (int,)) else list(dim2squish)

	out = data
	for axis in axes:
		out = BlockCat(out, axis, dim2concat)
	return out

class _ShapeWrapper:
	"""
	Lightweight wrapper to provide .shape and .ndim without creating a full ndarray.
	"""
	def __init__(self, shape: tuple[int, ...]):
		self.shape = shape
		self.ndim = len(shape)


def invMultiBlockCat(
	concatenated: np.ndarray,
	original: np.ndarray,
	dim2squish: int | list[int],
	dim2concat: int,
	was_squeezed: bool = False
) -> np.ndarray:
	"""
	Reverse of multiBlockCat: given the array after multiBlockCat(original,...),
	reconstruct the original array by applying invBlockCat in reverse order.

	Parameters
	----------
	concatenated : np.ndarray
		The result of multiBlockCat(original, dim2squish, dim2concat),
		possibly squeezed (if was_squeezed=True).
	original : np.ndarray
		The array that was fed into multiBlockCat.
	dim2squish : int or list of int
		Axis or axes along which original was split.
	dim2concat : int
		The axis along which concatenations were done.
	was_squeezed : bool, optional
		If True, indicates the concatenated array had singleton dims removed
		via squeeze() after multiBlockCat.  Those dims will be re-inserted.

	Returns
	-------
	np.ndarray
		The reconstructed original array.
	"""
	# Normalize axes list
	axes = [dim2squish] if isinstance(dim2squish, int) else list(dim2squish)

	# Compute intermediate shapes
	shapes: list[tuple[int, ...]] = [original.shape]
	for axis in axes:
		prev = shapes[-1]
		i = axis if axis >= 0 else len(prev) + axis
		j = dim2concat if dim2concat >= 0 else len(prev) + dim2concat
		n_slices = prev[i]

		new_shape = list(prev)
		# after BlockCat: axis i becomes size=1
		new_shape[i] = 1
		# and axis j is multiplied by n_slices (unless i==j)
		if i != j:
			new_shape[j] = prev[j] * n_slices
		shapes.append(tuple(new_shape))

	# If the array was squeezed after multiBlockCat, restore singleton dims
	if was_squeezed:
		expected = shapes[-1]
		concatenated = concatenated.reshape(expected)

	# Invert in reverse order
	out = concatenated
	for axis, orig_shape in zip(reversed(axes), reversed(shapes[:-1])):
		wrapper = _ShapeWrapper(orig_shape)
		out = invBlockCat(out, wrapper, axis, dim2concat)
	return out


class BlockTransformer:
	"""
	Encapsulates BlockCat and its inverse variants.

	Example:
		transformer = BlockCatTransformer(data, dim2squish=[0, 2], dim2concat=1)
		transformed = transformer.multi_block_cat()
		recovered = transformer.inv_multi_block_cat()
	"""
	def __init__(
		self,
		data: np.ndarray,
		dim2squish: int | list[int],
		dim2concat: int,
		squeeze: bool=False
	):
		# Store original data and parameters
		self.original = data
		self.dim2squish = dim2squish
		self.dim2concat = dim2concat
		self.transformed = None
		self.was_squeezed=squeeze

	def multi_block_cat(self) -> np.ndarray:
		"""
		Apply multiBlockCat (possibly multiple axes) and cache result.
		"""
		self.transformed = multiBlockCat(
			self.original,
			self.dim2squish,
			self.dim2concat
		)
		if self.was_squeezed:
			self.transformed=self.transformed.squeeze()
		return self.transformed

	def inv_multi_block_cat(self,array=None,org_array=None) -> np.ndarray:
		"""
		Apply invMultiBlockCat and return the reconstructed original.
		"""
		if array is None:
			array=self.transformed
		if org_array is None:
			org_array=self.original
		return invMultiBlockCat(
			array,
			org_array,
			dim2squish=self.dim2squish,
			dim2concat=self.dim2concat,
			was_squeezed=self.was_squeezed,
		)


def getColorMap(name:str,size:int): 
	cmap=plt.get_cmap(name)
	cmap=cmap(np.linspace(0,1,size))
	return cmap[:,:3]

def jumping_mean(data,dim,num_to_combine):
	# if num_to_combine is infinite, return mean of entire data
	if num_to_combine==np.inf:
		return data.mean(axis=dim,keepdims=True)

	data=np.asanyarray(data)
	stepsize=np.round(data.shape[dim]/num_to_combine).astype(int)
	mean_data=[]
	for start in range(stepsize):
		ind=np.arange(start,data.shape[dim],stepsize).astype(int)
		currdata=data.take(ind,axis=dim)
		mean_data.append(currdata.mean(axis=dim))
	mean_data=np.stack(mean_data,axis=dim)
	return mean_data
def force_monotonic_increase(x,min_x=None,max_x=None):
	"""
	Auto-detect the largest strictly-increasing contiguous segment in a 1D array
	and adjust the endpoints so the entire vector follows the median slope of
	that segment.

	Parameters
	----------
	x : array-like, shape (N,)
		Input vector, usually near-linear but corrupted at ends.
	
	Returns
	-------
	x_fixed : ndarray, shape (N,)
		A copy of x with its first/last values replaced so that the full
		sequence follows the detected constant slope.
	"""
	x = np.asarray(x, dtype=float)
	N = x.size
	if N < 3:
		return x.copy()

	# 1) find diffs and mask of strictly increasing steps
	diffs = np.diff(x)
	good = diffs > 0

	# 2) find longest run of True in 'good'
	max_len = 0
	max_start = 0
	curr_len = 0
	curr_start = None
	for i, val in enumerate(good):
		if val:
			if curr_start is None:
				curr_start = i
				curr_len = 1
			else:
				curr_len += 1
		else:
			if curr_start is not None and curr_len > max_len:
				max_len, max_start = curr_len, curr_start
			curr_start, curr_len = None, 0
	# check last run
	if curr_start is not None and curr_len > max_len:
		max_len, max_start = curr_len, curr_start

	# 3) determine region and slope
	if max_len < 1:
		region_start, region_end = 0, N - 1
		slope = np.median(diffs)
	else:
		region_start = max_start
		region_end = max_start + max_len
		slope = np.median(diffs[region_start:region_end])

	# 4) build ideal straight line anchored at region_start
	anchor_idx = region_start
	anchor_val = x[anchor_idx]
	idx = np.arange(N)
	ideal = anchor_val + slope * (idx - anchor_idx)

	# 5) apply fixes only outside the detected region
	x_fixed = x.copy()
	x_fixed[:region_start] = ideal[:region_start]
	x_fixed[region_end + 1:] = ideal[region_end + 1:]
	
	x_fixed=np.clip(x_fixed,min_x,max_x)

	return x_fixed

def save_as_svg(fig,path=None,name=None):
	if path is None:
		path=r"/path/to/data/tempfolder"
		if name is None:
			name=r"temp"
		name=name+".svg"
		path=os.path.join(path,name)
	else:
		path=path+".svg"
	fig.savefig(path)


def save_as_png(fig, path=None, name=None, dpi=200, **savefig_kwargs):
	"""
	Save a matplotlib Figure as PNG.

	Interface mirrors the provided save_as_svg:
	  - If path is None: use default folder and (optional) name.
	  - If path is given: treat it as a full file path *without* extension.

	Parameters
	----------
	fig : matplotlib.figure.Figure
		Figure to save.
	path : str | None
		If None, uses default folder. If provided, treated as full path w/o extension.
	name : str | None
		Filename (without extension) when path is None.
	dpi : int, default 200
		PNG resolution.
	**savefig_kwargs :
		Passed through to fig.savefig (e.g., bbox_inches='tight').

	Returns
	-------
	str
		The full path to the saved PNG.
	"""
	if path is None:
		path = r"/path/to/data/tempfolder"
		if name is None:
			name = r"temp"
		filename = name + ".png"
		os.makedirs(path, exist_ok=True)
		full_path = os.path.join(path, filename)
	else:
		# Treat 'path' as full file path WITHOUT extension (same as your SVG example)
		full_path = path + ".png"
		parent = os.path.dirname(full_path)
		if parent:
			os.makedirs(parent, exist_ok=True)

	fig.savefig(full_path, dpi=dpi, **savefig_kwargs)
	return full_path

	
	
class timer:
	def __init__(self):
		self.start = time.perf_counter()
	def __call__(self):
		elapsed = time.perf_counter() - self.start
		hours, rem = divmod(elapsed, 3600)
		minutes, seconds = divmod(rem, 60)
		parts = []
		if hours >= 1:
			parts.append(f"{int(hours)} h")
		if minutes >= 1:
			parts.append(f"{int(minutes)} min")
		# Always include seconds (with milliseconds)
		parts.append(f"{seconds:0.3f} sec")
		output = ", ".join(parts)
		print(output)
		

def plot_circular_colormap(cmap_array, clockwise=False, start_angle=0,
									  tick_interval=45, ax=None, radius=1.0, width=0.3):
	"""
	Plots a circular colormap as a ring of colored wedges and adds radial ticks.

	Parameters:
	- cmap_array: (n,3) array of RGB colors in [0,1].
	- clockwise: bool. If True, wedges advance clockwise; if False, counter-clockwise.
	- start_angle: starting angle in degrees (0° is along +x axis).
	- tick_interval: degrees between ticks (default 45°).
	- ax: optional Matplotlib Axes. If None, uses current axes (plt.gca()).
	- radius: outer radius of the color ring.
	- width: thickness of the ring.
	"""
	n = cmap_array.shape[0]
	segment = 360.0 / n
	if ax is None:
		ax = plt.gca()

	# Draw wedges
	for i, color in enumerate(cmap_array):
		if clockwise:
			theta1 = start_angle - (i + 1) * segment
			theta2 = start_angle - i * segment
		else:
			theta1 = start_angle + i * segment
			theta2 = start_angle + (i + 1) * segment

		wedge = Wedge((0, 0), radius, theta1, theta2,
					  width=width, facecolor=color, edgecolor='none')
		ax.add_patch(wedge)

	# Draw ticks & labels
	tick_length  = 0.1
	label_offset = 0.15
	for angle in np.arange(0, 360, tick_interval):
		tick_angle = start_angle - angle if clockwise else start_angle + angle
		rad = np.deg2rad(tick_angle)

		x0 = (radius - width/2) * np.cos(rad)
		y0 = (radius - width/2) * np.sin(rad)
		x1 = (radius + tick_length) * np.cos(rad)
		y1 = (radius + tick_length) * np.sin(rad)
		ax.plot([x0, x1], [y0, y1], color='k', lw=1)

		xt = (radius + label_offset) * np.cos(rad)
		yt = (radius + label_offset) * np.sin(rad)
		ax.text(xt, yt, f"{int(angle)}°", ha='center', va='center', fontsize=10)

	ax.set_aspect('equal')
	lim = radius + label_offset + tick_length
	ax.set_xlim(-lim, lim)
	ax.set_ylim(-lim, lim)
	ax.axis('off')
	return ax

def contig(array):
	return np.ascontiguousarray(array,dtype=np.float64)



def plot_mean_sem(A, x=None, ax=None, label=None, ci=1.0, **plot_kw):
	"""
	Plot row-wise mean ± ci*SEM for a 2D array A of shape (rows, cols),
	reducing across the SECOND dimension (axis=1, columns).

	Result: one point per row (length = A.shape[0]).
	- NaNs are ignored per row.
	- ci=1.0 -> mean ± 1*SEM; use 1.96 for ~95% CI.
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	A = np.asarray(A, dtype=float)
	if A.ndim != 2:
		raise ValueError("A must be 2D")

	rows = A.shape[0]

	# x must have length = number of rows
	if x is None:
		x = np.arange(rows)
	else:
		x = np.asarray(x).ravel()
		if x.size != rows:
			raise ValueError(f"x must have length {rows}, got {x.size}")

	if ax is None:
		ax = plt.gca()

	# counts, mean, std, sem per ROW (reducing across columns)
	n    = np.sum(~np.isnan(A), axis=1)
	mean = np.nanmean(A, axis=1)
	std  = np.nanstd(A, axis=1, ddof=1)

	sem = np.zeros_like(mean)
	valid = n > 1
	sem[valid] = std[valid] / np.sqrt(n[valid])
	# if n<=1, sem stays 0 (no band)

	(line,) = ax.plot(x, mean, label=label, **plot_kw)
	ax.fill_between(x, mean - ci*sem, mean + ci*sem,
					alpha=0.25, linewidth=0, color=line.get_color())

	# ax.set_xlabel("Row index" if x is None else ax.get_xlabel())
	# ax.set_ylabel("Value")
	# if label is not None:
	# 	ax.legend()
	return ax

####################### statistical test
def plot_statTest(
	A, x=None, ax=None, base_row=0,
	test='welch',                   # 'welch','ttest','mannwhitney','paired','wilcoxon'
	alpha=0.05, p_adjust=None, show_ns=False,
	offset_points=(6, 4), star_levels=(0.05, 0.01, 0.001),
	text_kw=None,
	alternative='two-sided',        # used for MWU and Wilcoxon
	zero_method='wilcox'            # how to handle zero diffs in Wilcoxon ('wilcox','pratt','zsplit')
):
	"""
	Compare each row to `base_row` and annotate significance at (x[i], mean(A[i,:])).
	No mean/SEM lines are drawn—only significance markers (*/**/*** and ↑/↓).

	test:
	  - 'welch'       : independent Welch t-test (unequal variances)
	  - 'ttest'       : independent Student t-test (equal variances)
	  - 'mannwhitney' : independent Mann–Whitney U (two-sided)
	  - 'paired'      : paired t-test on within-subject differences (row i - base_row)
	  - 'wilcoxon'    : Wilcoxon signed-rank on paired differences

	Returns
	-------
	ax, results dict with keys: 'means','p_raw','p_adj','effect_d','direction'
	  effect_d:
		- t-tests: Cohen's d (independent) or d_z (paired; mean(diff)/sd(diff))
		- MWU: rank-biserial proxy (independent)
		- Wilcoxon: rank-biserial correlation for signed-rank
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	try:
		from scipy import stats
	except Exception as e:
		raise ImportError("plot_statTest requires SciPy (scipy.stats).") from e

	A = np.asarray(A, dtype=float)
	if A.ndim != 2:
		raise ValueError("A must be 2D")
	rows, cols = A.shape

	if x is None:
		x = np.arange(rows)
	else:
		x = np.asarray(x).ravel()
		if x.size != rows:
			raise ValueError(f"x must have length {rows}, got {x.size}")

	if ax is None:
		ax = plt.gca()

	# Row means for placing annotations
	mean = np.nanmean(A, axis=1)

	# Baseline row vector (keep NaNs; we’ll pair per-column as needed)
	b0 = A[base_row, :]

	p_raw = np.full(rows, np.nan, dtype=float)
	effect_d = np.full(rows, np.nan, dtype=float)
	direction = np.zeros(rows, dtype=int)

	# Helpers
	def cohend_ind(x, y):
		x = x[~np.isnan(x)]; y = y[~np.isnan(y)]
		nx, ny = x.size, y.size
		if nx < 2 or ny < 2: return np.nan
		sx = np.nanstd(x, ddof=1); sy = np.nanstd(y, ddof=1)
		df = (nx - 1) + (ny - 1)
		if df <= 0: return np.nan
		sp2 = ((nx - 1)*sx*sx + (ny - 1)*sy*sy) / df
		if sp2 <= 0 or not np.isfinite(sp2): return np.nan
		return (np.nanmean(x) - np.nanmean(y)) / np.sqrt(sp2)

	def cohend_paired(d):
		d = d[~np.isnan(d)]
		if d.size < 2: return np.nan
		sd = np.nanstd(d, ddof=1)
		if sd == 0 or not np.isfinite(sd): return np.nan
		return np.nanmean(d) / sd

	for i in range(rows):
		if i == base_row:
			continue

		ai = A[i, :]

		if test in ('paired', 'wilcoxon'):
			# Paired design: use within-column differences
			mask = (~np.isnan(ai)) & (~np.isnan(b0))
			d = ai[mask] - b0[mask]
			if d.size == 0:
				continue

			mdiff = np.nanmean(d)
			direction[i] = 0 if (not np.isfinite(mdiff) or np.isclose(mdiff, 0.0)) else (1 if mdiff > 0 else -1)

			if test == 'paired':
				if d.size < 2:
					continue
				# t-test on differences vs 0
				res = stats.ttest_1samp(d, 0.0, nan_policy='omit')
				p_raw[i] = res.pvalue
				effect_d[i] = cohend_paired(d)

			else:  # 'wilcoxon'
				# SciPy wilcoxon handles zeros per zero_method; get W+ to compute rank-biserial
				try:
					res = stats.wilcoxon(d, alternative=alternative, zero_method=zero_method, method='auto')
				except TypeError:
					# Older SciPy fallback (no 'method' arg)
					res = stats.wilcoxon(d, alternative=alternative, zero_method=zero_method)
				p_raw[i] = res.pvalue

				# Rank-biserial correlation for signed-rank:
				# W+ is the returned statistic; n_eff = number of nonzero diffs if zero_method='wilcox'
				if zero_method == 'wilcox':
					n_eff = np.count_nonzero(d)
				elif zero_method == 'pratt':
					# Zeros included in ranking; SciPy drops zero ranks from sum; use nonzero count for T
					n_eff = np.count_nonzero(d)
				else:  # 'zsplit' splits zeros; treat them as half-positive/half-negative -> still use nonzero count
					n_eff = np.count_nonzero(d)

				if n_eff >= 1:
					T = n_eff * (n_eff + 1) / 2.0
					Wplus = float(res.statistic)
					r_rb = 2.0 * (Wplus / T) - 1.0
					effect_d[i] = r_rb
				else:
					effect_d[i] = np.nan

		else:
			# Independent-samples branch
			xg = ai[~np.isnan(ai)]
			yg = b0[~np.isnan(b0)]
			if xg.size == 0 or yg.size == 0:
				continue

			mdiff = np.nanmean(xg) - np.nanmean(yg)
			direction[i] = 0 if (not np.isfinite(mdiff) or np.isclose(mdiff, 0.0)) else (1 if mdiff > 0 else -1)

			if test in ('welch', 'ttest'):
				eq = (test == 'ttest')
				if xg.size < 2 or yg.size < 2:
					continue
				res = stats.ttest_ind(xg, yg, equal_var=eq, nan_policy='omit')
				p_raw[i] = res.pvalue
				effect_d[i] = cohend_ind(xg, yg)
			elif test == 'mannwhitney':
				try:
					res = stats.mannwhitneyu(xg, yg, alternative=alternative, method='auto')
				except TypeError:
					res = stats.mannwhitneyu(xg, yg, alternative=alternative)
				p_raw[i] = res.pvalue
				# rank-biserial for independent samples
				nx, ny = xg.size, yg.size
				U = float(res.statistic)
				effect_d[i] = 1.0 - (2.0 * U) / (nx * ny)
			else:
				raise ValueError("test must be one of {'welch','ttest','mannwhitney','paired','wilcoxon'}")

	# Multiple-comparisons correction across tested rows
	def _bonf(p):
		q = p.copy()
		m = np.sum(~np.isnan(q))
		if m > 0:
			q = np.minimum(q * m, 1.0)
		return q

	def _fdr_bh(p):
		p = p.copy()
		idx = np.where(~np.isnan(p))[0]
		if idx.size == 0:
			return p
		ps = p[idx]
		order = np.argsort(ps)
		ps_sorted = ps[order]
		m = ps_sorted.size
		ranks = np.arange(1, m + 1)
		q_sorted = ps_sorted * m / ranks
		for k in range(m - 2, -1, -1):
			q_sorted[k] = min(q_sorted[k], q_sorted[k + 1])
		q_sorted = np.minimum(q_sorted, 1.0)
		p_adj = p.copy()
		p_adj[idx[order]] = q_sorted
		return p_adj

	if p_adjust is None:
		p_adj = p_raw.copy()
	elif str(p_adjust).lower() == 'bonferroni':
		p_adj = _bonf(p_raw)
	elif str(p_adjust).lower() in ('fdr', 'fdr_bh', 'bh'):
		p_adj = _fdr_bh(p_raw)
	else:
		raise ValueError("p_adjust must be None, 'bonferroni', or 'fdr_bh'")

	# Stars
	def p_to_stars(p):
		if not np.isfinite(p): return ''
		if p <= star_levels[2]: return '***'
		if p <= star_levels[1]: return '**'
		if p <= star_levels[0]: return '*'
		return ''

	if text_kw is None:
		text_kw = {}
	dx, dy = offset_points

	for i in range(rows):
		if i == base_row:
			continue
		yi = mean[i]
		if not np.isfinite(yi):
			continue
		stars = p_to_stars(p_adj[i])
		if stars or show_ns:
			arrow = ''
			if stars or show_ns:
				arrow = '↑' if direction[i] > 0 else ('↓' if direction[i] < 0 else '')
			txt = stars + (arrow if stars else ('ns' if show_ns else ''))
			ax.annotate(txt, (x[i], yi),
						textcoords='offset points', xytext=(dx, dy),
						ha='left', va='bottom', **text_kw)

	results = {
		'means': mean,
		'p_raw': p_raw,
		'p_adj': p_adj,
		'effect_d': effect_d,
		'direction': direction
	}
	return ax, results

######################################


def set_axis_limits(ax, X, margin=0.1):
	for i, setter in enumerate([ax.set_xlim, ax.set_ylim]):
		dmin, dmax = X[:, i].min(), X[:, i].max()
		m = margin * (dmax - dmin)
		setter(dmin - m, dmax + m)



from pathlib import Path
from matplotlib import pyplot as plt
import shutil
import re
class FigSaver:
	"""
	Save matplotlib figures sequentially to a folder, then stitch them into a GIF.

	Parameters
	----------
	folder : str | Path
		Destination folder (created if it doesn't exist).
	prefix : str, default 'fig'
		Filename prefix, e.g., fig_0001.png
	start : int, default 1
		Starting index for sequential filenames.
	zpad : int, default 4
		Zero-padding width for the index.
	ext : str, default 'png'
		Image file extension (without dot).
	dpi : int, default 200
		Resolution for saved images.
	transparent : bool, default False
		Save with transparent background.
	bbox_inches : 'tight' | None, default 'tight'
		Passed to savefig to trim margins.
	pad_inches : float, default 0.02
		Padding (inches) when using bbox_inches='tight'.

	Methods
	-------
	save(fig=None, name=None, close=False, **kwargs) -> Path
		Save a figure (or current figure if None). Returns saved Path.
	save_all(figs, close=False, **kwargs) -> list[Path]
		Save a list/iterable of figures.
	to_gif(gif_path=None, fps=10, loop=0, clean=False, use_saved_order=True) -> Path
		Make a GIF from saved frames. Optionally delete the image folder afterward.
	"""

	def __init__(self, folder, prefix='fig', start=1, zpad=4,
				 ext='png', dpi=200, transparent=False,
				 bbox_inches='tight', pad_inches=0.02):
		self.folder = Path(folder)
		self.folder.mkdir(parents=True, exist_ok=True)

		# Prefer imageio; fall back to PIL if needed
		try:
			import imageio.v3 as iio  # modern API
			_USE_IMAGEIO = True
		except Exception:
			from PIL import Image
			_USE_IMAGEIO = False

		self.prefix = prefix
		self.counter = int(start)
		self.zpad = int(zpad)
		self.ext = ext.lstrip('.')

		self._default_savefig = dict(
			dpi=dpi,
			transparent=transparent,
			bbox_inches=bbox_inches,
			pad_inches=pad_inches
		)

		# Keep track of what we saved (order-preserving)
		self._saved_paths = []

	def _next_filename(self):
		name = f"{self.prefix}_{self.counter:0{self.zpad}d}.{self.ext}"
		self.counter += 1
		return self.folder / name

	def save(self, fig=None, name=None, close=False, **kwargs):
		fig = plt.gcf() if fig is None else fig
		path = (self.folder / name) if name else self._next_filename()
		if path.suffix.lower() != f".{self.ext}":
			path = path.with_suffix(f".{self.ext}")

		savefig_kwargs = {**self._default_savefig, **kwargs}
		fig.savefig(path, **savefig_kwargs)
		if close:
			plt.close(fig)

		self._saved_paths.append(path)
		return path

	def save_all(self, figs, close=False, **kwargs):
		return [self.save(f, close=close, **kwargs) for f in figs]

	# -------- GIF creation --------
	@staticmethod
	def _natural_key(p: Path):
		"""Sort key that handles numbers in filenames naturally (e.g., 2 < 10)."""
		parts = re.split(r'(\d+)', p.name)
		return [int(s) if s.isdigit() else s.lower() for s in parts]

	def _collect_frames(self, use_saved_order: bool):
		if use_saved_order and self._saved_paths:
			frames = [p for p in self._saved_paths if p.exists()]
		else:
			# Fall back to discovering images in the folder (natural sort)
			frames = sorted(self.folder.glob(f"*.{self.ext}"), key=self._natural_key)
		return frames

	def to_gif(self, gif_path=None, fps=10, loop=0, clean=False, use_saved_order=True):
		"""
		Create an animated GIF from the saved images.

		Parameters
		----------
		gif_path : str | Path | None
			Full path (including filename) to the output GIF.
			If None, writes '<parent>/<foldername>.gif'.
		fps : int, default 10
			Frames per second (speed of animation).
		loop : int, default 0
			Number of times the GIF will loop. 0 = infinite looping.
		clean : bool, default False
			If True, delete the entire image folder AFTER writing the GIF.
		use_saved_order : bool, default True
			If True, use the order images were saved in this FigSaver session.
			If False (or if history is empty), discover and sort files by name.

		Returns
		-------
		Path
			Path to the created GIF.
		"""
		frames = self._collect_frames(use_saved_order=use_saved_order)
		if not frames:
			raise RuntimeError("No image frames found to build a GIF.")

		# Output path
		if gif_path is None:
			gif_path = self.folder.parent / f"{self.folder.name}.gif"
		gif_path = Path(gif_path)
		if gif_path.suffix.lower() != ".gif":
			gif_path = gif_path.with_suffix(".gif")
		gif_path.parent.mkdir(parents=True, exist_ok=True)

		# Safety: don't allow deleting the folder if the GIF is placed inside it
		# (otherwise you'd delete your result)
		if clean:
			try:
				folder_resolved = self.folder.resolve()
				gif_parent_resolved = gif_path.resolve().parent
				if folder_resolved == gif_parent_resolved or str(gif_path.resolve()).startswith(str(folder_resolved) + str(Path.sep)):
					raise ValueError(
						"gif_path points inside the image folder while clean=True. "
						"Write the GIF outside the folder if you want to delete the folder afterward."
					)
			except Exception:
				# If resolve() fails for any reason, fall back to a simple string check
				if str(gif_path).startswith(str(self.folder)):
					raise ValueError(
						"gif_path points inside the image folder while clean=True."
					)

		# Build GIF
		if _USE_IMAGEIO:
			# Read all frames as arrays; duration is seconds per frame
			imgs = [iio.imread(p) for p in frames]
			iio.imwrite(gif_path, imgs, plugin="pillow", duration=1.0 / fps, loop=loop)
		else:
			# PIL fallback
			images = [Image.open(p).convert("RGBA") for p in frames]
			images[0].save(
				gif_path,
				save_all=True,
				append_images=images[1:],
				duration=int(1000 / fps),  # ms per frame
				loop=loop,
				disposal=2,
			)

		# Optionally delete the image folder
		if clean:
			shutil.rmtree(self.folder)

		return gif_path



def polygon_area(X, *, signed: bool = False) -> float:
	"""
	Compute the area of a simple 2D polygon from its vertices.

	Parameters
	----------
	X : array_like, shape (n, 2)
		Polygon vertices in order (clockwise or counterclockwise).
		The polygon may be open (first point != last point); it will
		be closed internally.
	signed : bool, optional (default: False)
		If True, return signed area (CCW positive, CW negative).
		If False, return absolute area.

	Returns
	-------
	float
		Polygon area.
	"""
	P = np.asarray(X, dtype=float)
	if P.ndim != 2 or P.shape[1] != 2:
		raise ValueError("X must be an (n, 2) array of 2D points.")
	if len(P) < 3:
		return 0.0

	# Ensure closed ring
	if not np.allclose(P[0], P[-1]):
		P = np.vstack([P, P[0]])

	# Shoelace sum
	cross = P[:-1, 0] * P[1:, 1] - P[1:, 0] * P[:-1, 1]
	area = 0.5 * cross.sum()
	return area if signed else abs(area)


# Optional: helper to enforce CCW ordering
def ensure_ccw(X):
	"""Return a copy of X ordered counterclockwise."""
	return np.array(X)[::-1].copy() if polygon_area(X, signed=True) < 0 else np.array(X)



import numpy as np

def polygon_area_union(X) -> float:
	from shapely.geometry import Polygon
	P = np.asarray(X, float)
	if len(P) < 3: 
		return 0.0
	if not np.allclose(P[0], P[-1]):
		P = np.vstack([P, P[0]])
	poly = Polygon(P)          # may be invalid (self-intersecting)
	fixed = poly.buffer(0)     # resolves self-intersections into valid geometry
	return fixed.area


def angle_color(angles, cmap_name='hsv',nStep=256):
	"""
	Map angles in radians to RGB colors using a circular colormap.

	Parameters
	----------
	angles : array-like
		Input angles in radians.
	cmap_name : str, optional
		Name of a Matplotlib colormap (default: 'hsv').
	nStep : int, optional
		Number of discrete colors to sample from the colormap (default: 256).

	Returns
	-------
	rgb_colors : ndarray, shape (..., 3)
		RGB colors corresponding to input angles.
	"""
	import numpy as np
	from matplotlib import cm

	angles = np.asarray(angles)
	cmap = cm.get_cmap(cmap_name, nStep)
	normed = (angles % (2 * np.pi)) / (2 * np.pi)  # Normalize to [0, 1)
	indices = np.round(normed * (nStep - 1)).astype(int)
	rgb_colors = cmap(indices)[..., :3]  # Discard alpha channel if present
	return rgb_colors



def fix_strictly_increasing(X, axis=-1, eps=None):
	"""
	Make X strictly increasing along `axis` by anchoring the longest strictly-increasing
	contiguous segment and minimally adjusting only the values outside that segment.

	Parameters
	----------
	X : np.ndarray
		Input array. Along `axis`, values are intended to be strictly increasing,
		but edges may violate monotonicity.
	axis : int, default -1
		Axis along which to enforce strict increase.
	eps : float or None, default None
		Minimum step. If None, uses a small dtype-dependent epsilon (float).
		For integer dtypes, X is promoted to float64 and eps is set automatically.

	Returns
	-------
	Y : np.ndarray
		Copy of X with strictly increasing values along `axis`. The longest
		strictly-increasing segment found is preserved exactly; only the parts
		outside that segment are adjusted, with minimal changes to satisfy monotonicity.
	anchor_slices : tuple
		A tuple of slice objects marking the anchored segment along `axis`
		(useful for inspection/plotting).
	"""
	A = np.asarray(X)
	# Ensure we work in float with a sensible epsilon
	if A.dtype.kind in ("i", "u"):
		work = A.astype(np.float64, copy=True)
		local_eps = 1e-12 if eps is None else float(eps)
	else:
		work = A.astype(float, copy=True)
		if eps is None:
			finfo = np.finfo(work.dtype)
			# a bit larger than machine epsilon to avoid numerical ties
			local_eps = 10 * finfo.eps
		else:
			local_eps = float(eps)

	# Move target axis to last, then flatten the rest to rows
	work = np.moveaxis(work, axis, -1)
	leading_shape, L = work.shape[:-1], work.shape[-1]
	work2d = work.reshape(-1, L)  # rows of length L

	def longest_increasing_segment(v):
		"""Return (s, e) inclusive indices of longest strictly increasing contiguous segment."""
		d = np.diff(v)
		inc = d > 0
		if inc.size == 0:
			# length-1 vector: anchor the only element
			return 0, 0

		z = np.r_[False, inc, False]
		changes = np.diff(z.astype(int))
		starts = np.flatnonzero(changes == 1)    # indices in inc where runs start
		ends   = np.flatnonzero(changes == -1)   # indices in inc where runs end
		if starts.size == 0:
			# No strictly increasing step; anchor a single element with largest value
			j = int(np.argmax(v))
			return j, j

		lengths = ends - starts                  # run lengths in 'inc' (steps)
		# break ties by total rise over the run (prefer bigger overall increase)
		rises = np.array([d[s:s+ln].sum() for s, ln in zip(starts, lengths)])
		best = np.lexsort((-rises, -lengths))[0]  # maximize length, then rise
		s = starts[best]
		e = s + lengths[best]                    # inclusive end in v
		return int(s), int(e)

	def repair_with_anchor(v, s, e, eps_):
		"""Fix values outside [s,e] so the whole vector is strictly increasing; keep [s,e] untouched."""
		out = v.copy()
		n = out.size
		# Left side: walk downwards, cap to just below next
		for i in range(s - 1, -1, -1):
			allowed = out[i + 1] - eps_
			if out[i] >= allowed:
				out[i] = allowed
		# Right side: walk upwards, lift to just above previous
		for i in range(e + 1, n):
			allowed = out[i - 1] + eps_
			if out[i] <= allowed:
				out[i] = allowed
		return out

	anchors = []
	for r in range(work2d.shape[0]):
		v = work2d[r]
		s, e = longest_increasing_segment(v)
		work2d[r] = repair_with_anchor(v, s, e, local_eps)
		anchors.append((s, e))

	# Restore original shape/axis
	out = work2d.reshape(*leading_shape, L)
	out = np.moveaxis(out, -1, axis)

	# Build a generic anchor_slices describing the anchored span along `axis`
	# (Note: for multidimensional inputs, each row can have a different anchor.)
	# Here we return per-row anchors (list of tuples), since that’s the most informative.
	return out


###########################################################################
# Statistical tests for monotonic trends in repeated-measures data
def _per_mouse_slopes(A, x):
	"""
	Compute per-mouse OLS slopes of score ~ x.
	A: (k_levels, n_mice), x: (k_levels,)
	Returns slopes (n_mice,), mask of mice actually used (>=2 valid points).
	"""
	x = np.asarray(x, float).reshape(-1)
	A = np.asarray(A, float)
	k, n = A.shape
	slopes = np.full(n, np.nan)
	for j in range(n):
		y = A[:, j]
		msk = np.isfinite(x) & np.isfinite(y)
		if msk.sum() >= 2:
			# simple OLS slope with intercept
			X = np.column_stack([np.ones(msk.sum()), x[msk]])
			beta, *_ = np.linalg.lstsq(X, y[msk], rcond=None)
			slopes[j] = beta[1]
	used = np.isfinite(slopes)
	return slopes[used], used

def _pages_L_test(A_sorted):
	"""
	Page's L trend test for repeated measures with ordered levels (rows).
	A_sorted: (k, n) with rows already sorted by ascending x.
	Returns dict with L, z, p_increasing, p_decreasing, p_two_sided.
	Drops mice (columns) that contain any NaN.
	"""
	A_sorted = np.asarray(A_sorted, float)
	k, n = A_sorted.shape
	# keep complete-case mice (no NaN across levels)
	ok = np.all(np.isfinite(A_sorted), axis=0)
	Y = A_sorted[:, ok]
	k, n2 = Y.shape
	if n2 < 2 or k < 2:
		return {"L": np.nan, "z": np.nan, "p_increasing": np.nan,
				"p_decreasing": np.nan, "p_two_sided": np.nan, "n_mice": n2}
	# rank within each mouse (column); smaller = rank 1
	# average ranks for ties
	R = np.zeros_like(Y)
	for j in range(n2):
		col = Y[:, j]
		# argsort twice trick for ranks with ties -> use scipy if you prefer
		order = np.argsort(col, kind="mergesort")
		ranks = np.empty(k, float)
		ranks[order] = np.arange(1, k+1)
		# tie handling (average ranks)
		# find ties in sorted order
		s = col[order]
		start = 0
		while start < k:
			end = start + 1
			while end < k and s[end] == s[start]:
				end += 1
			if end - start > 1:
				avg = (start + 1 + end) / 2.0
				ranks[order[start:end]] = avg
			start = end
		R[:, j] = ranks

	Rj = R.sum(axis=1)                 # sum of ranks for each level j
	j_idx = np.arange(1, k+1, dtype=float)
	L = float(np.dot(j_idx, Rj))       # Page's L statistic

	# Null moments (no ordered effect), large-sample normal approx
	# E[L] = n * k * (k+1)^2 / 4
	# Var[L] = n * k^2 * (k+1)^2 * (k-1) / 144
	E = n2 * k * (k + 1) ** 2 / 4.0
	Var = n2 * (k ** 2) * ((k + 1) ** 2) * (k - 1) / 144.0
	z = (L - E) / np.sqrt(Var)

	# One-sided p for INCREASING trend is 1 - Phi(z)
	# For DECREASING trend use Phi(z)
	from math import erf, sqrt
	def Phi(u):  # CDF of N(0,1)
		return 0.5 * (1 + erf(u / sqrt(2)))
	p_inc = 1 - Phi(z)
	p_dec = Phi(z)
	p_two = 2 * min(p_inc, p_dec)
	return {"L": L, "z": z, "p_increasing": p_inc, "p_decreasing": p_dec,
			"p_two_sided": p_two, "n_mice": n2}

def slope_statistics(A, x, alternative="two-sided"):
	"""
	Test whether the score changes monotonically.
	- Mixed-effects (parametric): y ~ x + (1|mouse)
	- Page's L (nonparametric repeated-measures monotonic trend)
	- Per-mouse slopes -> Wilcoxon signed-rank vs 0 (and one-sample t-test)

	Parameters
	----------
	A : array (k_levels, n_mice)
	x : array (k_levels,)
	alternative : 'two-sided', 'increasing', or 'decreasing' for the slope tests

	Returns: dict of results
	"""
	import numpy as np
	x = np.asarray(x, float).reshape(-1)
	A = np.asarray(A, float)
	assert A.shape[0] == x.size, "A rows must match len(x)"

	# sort rows by ascending x (safer if x not pre-sorted)
	idx = np.argsort(x)
	x_sorted = x[idx]
	A_sorted = A[idx, :]

	# 1) Per-mouse slopes
	slopes, used = _per_mouse_slopes(A_sorted, x_sorted)
	n_mice = slopes.size

	from scipy.stats import ttest_1samp, wilcoxon, binomtest
	# Decide tail for slope-based tests
	if alternative not in ("two-sided", "increasing", "decreasing"):
		raise ValueError("alternative must be 'two-sided' | 'increasing' | 'decreasing'")
	# One-sample t-test on slopes
	tt = ttest_1samp(slopes, popmean=0.0, alternative=alternative if alternative!="two-sided" else "two-sided")
	# Wilcoxon signed-rank (requires at least one nonzero diff)
	try:
		wx = wilcoxon(slopes, zero_method="wilcox", alternative=alternative if alternative!="two-sided" else "two-sided", mode="auto")
	except Exception:
		wx = None
	# Simple sign test (very robust)
	n_pos = np.sum(slopes > 0)
	n_neg = np.sum(slopes < 0)
	n_use = n_pos + n_neg
	sign_p = binomtest(n_pos, n_use, p=0.5, alternative="two-sided" if alternative=="two-sided"
					   else ("greater" if alternative=="increasing" else "less")).pvalue if n_use>0 else np.nan

	# 2) Page's L test (monotone trend across ordered levels; repeated measures)
	page = _pages_L_test(A_sorted)

	# 3) Mixed-effects model with robust fallbacks
	mixed = None
	try:
		import numpy as np
		import pandas as pd
		import statsmodels.formula.api as smf
		import statsmodels.api as sm

		k, n = A_sorted.shape
		ys = A_sorted.reshape(-1)
		xs = np.repeat(x_sorted, n)
		mice = np.tile(np.arange(n), k)
		mask = np.isfinite(ys) & np.isfinite(xs)
		df = pd.DataFrame({"y": ys[mask], "x": xs[mask], "mouse": mice[mask]})

		# ----- Pre-filter: at least 2 usable points AND nonzero y variance within mouse -----
		def _mouse_ok(d):
			return (d["x"].nunique() >= 2) and (np.nanstd(d["y"]) > 0)
		df = df.groupby("mouse", group_keys=False).filter(_mouse_ok)

		# If still too few data, skip
		if df["mouse"].nunique() >= 2 and len(df) >= 6:
			# ----- Rescale x to [0,1] for stability & interpretability -----
			x_min, x_max = df["x"].min(), df["x"].max()
			if np.isfinite(x_min) and np.isfinite(x_max) and x_max > x_min:
				df = df.copy()
				df["x_s"] = (df["x"] - x_min) / (x_max - x_min)
				xvar = "x_s"
				x_scaled_01 = True
			else:
				xvar = "x"
				x_scaled_01 = False

			# ----- Try MixedLM with a few optimizer / REML combos -----
			last_err = None
			got = None
			for reml in (True, False):
				for method in ("lbfgs", "nm", "cg"):
					try:
						model = smf.mixedlm(f"y ~ {xvar}", df, groups=df["mouse"])
						fit = model.fit(method=method, reml=reml, disp=False)
						if np.isfinite(float(fit.params.get(xvar, np.nan))):
							beta = float(fit.params[xvar])
							pval = float(fit.pvalues[xvar])
							ci = fit.conf_int().loc[xvar].to_numpy().astype(float).tolist()
							mixed = {
								"beta_x": beta,
								"p": pval,
								"ci95": ci,
								"n_mice": int(df["mouse"].nunique()),
								"n_obs": int(len(df)),
								"x_scaled_01": x_scaled_01,
								"estimator": "MixedLM",
								"optimizer": method,
								"reml": reml,
							}
							got = True
							break
					except Exception as e:
						last_err = str(e)
				if got:
					break

			# ----- Fallback: Fixed-effects OLS with cluster-robust SE -----
			if mixed is None:
				try:
					fe = smf.ols(f"y ~ {xvar} + C(mouse)", data=df).fit(
						cov_type="cluster", cov_kwds={"groups": df["mouse"]}
					)
					beta = float(fe.params[xvar])
					se = float(fe.bse[xvar])
					pval = float(fe.pvalues[xvar])
					# Normal-based 95% CI (cluster-robust SE)
					from scipy.stats import norm
					z = norm.ppf(0.975)
					ci = [beta - z * se, beta + z * se]
					mixed = {
						"beta_x": beta,
						"p": pval,
						"ci95": ci,
						"n_mice": int(df["mouse"].nunique()),
						"n_obs": int(len(df)),
						"x_scaled_01": x_scaled_01,
						"estimator": "FE_cluster",
					}
				except Exception as e:
					print(f"MixedLM and FE fallback failed: {e}")
					mixed = None
		else:
			mixed = None

	except Exception as e:
		print(f"Error caught (skipping mixed): {e}")
		mixed = None


	# Direction guess from median slope
	med_slope = np.nanmedian(slopes) if slopes.size else np.nan
	direction = "increasing" if med_slope > 0 else ("decreasing" if med_slope < 0 else "flat")

	return {
		"x_sorted": x_sorted,
		"n_mice": int(n_mice),
		"slopes": slopes,
		"slope_summary": {
			"median": float(np.nanmedian(slopes)) if slopes.size else np.nan,
			"mean": float(np.nanmean(slopes)) if slopes.size else np.nan,
			"std": float(np.nanstd(slopes, ddof=1)) if slopes.size>1 else np.nan
		},
		"wilcoxon": ({"stat": float(wx.statistic), "p": float(wx.pvalue)} if wx is not None else None),
		"t_one_sample": {"stat": float(tt.statistic), "p": float(tt.pvalue)},
		"sign_test_p": float(sign_p),
		"page_L": page,     # includes one-sided and two-sided p-values
		"mixed_effects": mixed, # may be None if statsmodels unavailable
		"direction_estimate": direction,
	}

# --- OPTIONAL: tiny helper to format a short summary string ---
def summarize_trend_result(
	res,
	newline=False,
	include_wilcoxon=True,
	include_ttest=True,
	include_pageL=True,
	include_mixed=True,
):
	"""
	Build a summary string from analyze_wzvert_trend(res).

	Parameters
	----------
	res : dict
		Output from analyze_wzvert_trend(...)
	newline : bool, default False
		If True, each part is placed on its own line (joined by '\n').
		If False, parts are joined with ' | '.
	include_wilcoxon : bool
		Include Wilcoxon signed-rank test result.
	include_ttest : bool
		Include one-sample t-test result.
	include_pageL : bool
		Include Page's L result.
	include_mixed : bool
		Include mixed-effects model result (if available).

	Returns
	-------
	s : str
		Summary string.
	"""
	import numpy as np

	dir_est = res["direction_estimate"]
	s_med = res["slope_summary"]["median"]

	parts = [f"Direction: {dir_est} (median slope={s_med:.3g})"]

	if include_wilcoxon and res["wilcoxon"] is not None:
		parts.append(f"Wilcoxon p={res['wilcoxon']['p']:.3g}")

	if include_ttest and "t_one_sample" in res:
		parts.append(f"t-test p={res['t_one_sample']['p']:.3g}")

	if include_pageL and "page_L" in res:
		parts.append(f"Page's L two-sided p={res['page_L']['p_two_sided']:.3g}")

	if include_mixed and res.get("mixed_effects") is not None:
		mix = res["mixed_effects"]
		parts.append(
			f"MixedLM β_x={mix['beta_x']:.3g}, p={mix['p']:.3g}, "
			f"95% CI=({mix['ci95'][0]:.3g}, {mix['ci95'][1]:.3g})"
		)

	sep = "\n" if newline else " | "
	return sep.join(parts)



# # x: shape (k,), continuous w_zvert values
# # A: shape (k, n_mice), columns = mice
# res = slope_statistics(A, x, alternative="two-sided")  # or "increasing"/"decreasing" if you have a priori direction
# print(summarize_trend_result(res))
###########################################################################


###########################################################################
# statistics of 
# --- drop-in replacement: compare_to_baseline ---
import numpy as np
from scipy import stats

def _bh_fdr(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    idx = np.flatnonzero(np.isfinite(p))
    if len(idx) == 0: return out
    ps = p[idx]; m = len(ps)
    order = np.argsort(ps)
    ranks = np.empty_like(order); ranks[order] = np.arange(1, m+1)
    adj = ps * m / ranks
    adj_sorted = np.minimum.accumulate(adj[order[::-1]])[::-1]
    out[idx] = np.clip(adj_sorted[ranks-1], 0, 1)
    return out

def _holm(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    idx = np.flatnonzero(np.isfinite(p))
    if len(idx) == 0: return out
    ps = p[idx]; m = len(ps)
    order = np.argsort(ps)
    adj = np.empty_like(ps)
    for j, k in enumerate(order):
        adj[k] = ps[k] * (m - j)
    adj[order] = np.maximum.accumulate(adj[order])
    out[idx] = np.clip(adj, 0, 1)
    return out

def _bootstrap_ci(d, n_boot=5000, alpha=0.05, stat=np.nanmedian, random_state=None):
    rng = np.random.default_rng(random_state)
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    if d.size == 0: return (np.nan, np.nan)
    idx = rng.integers(0, d.size, size=(n_boot, d.size))
    boots = stat(d[idx], axis=1)
    lo = np.percentile(boots, 100*alpha/2)
    hi = np.percentile(boots, 100*(1 - alpha/2))
    return (lo, hi)

def _rank_biserial(diff):
    d = np.asarray(diff, float)
    d = d[np.isfinite(d) & (d != 0)]
    if d.size == 0: return np.nan
    ranks = stats.rankdata(np.abs(d), method="average")
    r_plus  = np.sum(ranks[d > 0])
    r_minus = np.sum(ranks[d < 0])
    return (r_plus - r_minus) / (r_plus + r_minus)

def compare_to_baseline(
	A, x,
	baseline_index=0,
	test="wilcoxon",            # 'wilcoxon' | 'ttest' | 'sign'
	alternative="two-sided",    # for wilcoxon/ttest
	p_adjust="fdr_bh",          # 'none' | 'holm' | 'fdr_bh'
	min_paired=3,
	n_boot=5000,
	random_state=None,
	store_matrix=False          # NEW: if True, include the whole A in the return
):
	"""
	A: (k_levels, n_mice). x: (k_levels,).
	Returns a dict with 'table' AND now also 'x_all', 'A_row_mean', 'A_row_median', etc.
	"""
	import pandas as pd

	A = np.asarray(A, float)
	x = np.asarray(x)
	k, n = A.shape
	if x.shape[0] != k: raise ValueError("x must have same length as A.shape[0].")
	if not (0 <= baseline_index < k): raise ValueError("baseline_index out of range.")

	rows = []
	base = A[baseline_index, :]

	for j in range(k):
		if j == baseline_index: continue
		a = A[j, :]
		mask = np.isfinite(base) & np.isfinite(a)
		nb = int(np.sum(mask))
		if nb < min_paired:
			rows.append({
				"row": j, "x": x[j], "n_paired": nb, "test": test,
				"stat": np.nan, "p": np.nan, "p_adj": np.nan,
				"median_diff": np.nan, "ci_low": np.nan, "ci_high": np.nan,
				"rank_biserial": np.nan, "direction": "insufficient data"
			})
			continue

		d = a[mask] - base[mask]
		med = np.nanmedian(d)
		ci_lo, ci_hi = _bootstrap_ci(d, n_boot=n_boot, stat=np.nanmedian, random_state=random_state)

		if test == "wilcoxon":
			res = stats.wilcoxon(a[mask], base[mask], alternative=alternative, zero_method='wilcox', correction=False)
			stat_val, pval = res.statistic, res.pvalue
		elif test == "ttest":
			res = stats.ttest_rel(a[mask], base[mask], alternative=alternative, nan_policy='omit')
			stat_val, pval = res.statistic, res.pvalue
		elif test == "sign":
			nonzero = d != 0
			n_eff = int(np.sum(nonzero))
			if n_eff == 0:
				stat_val, pval = np.nan, 1.0
			else:
				n_pos = int(np.sum(d[nonzero] > 0))
				pval = stats.binomtest(n_pos, n_eff, 0.5, alternative=("two-sided" if alternative=="two-sided" else alternative)).pvalue
				stat_val = n_pos
		else:
			raise ValueError("test must be 'wilcoxon', 'ttest', or 'sign'.")

		r_rb = _rank_biserial(d)
		rows.append({
			"row": j, "x": x[j], "n_paired": nb, "test": test,
			"stat": float(stat_val) if np.isfinite(stat_val) else np.nan,
			"p": float(pval) if np.isfinite(pval) else np.nan,
			"p_adj": np.nan,
			"median_diff": float(med), "ci_low": float(ci_lo), "ci_high": float(ci_hi),
			"rank_biserial": float(r_rb) if np.isfinite(r_rb) else np.nan,
			"direction": None
		})

	df = pd.DataFrame(rows).sort_values("x").reset_index(drop=True)

	# adjust p's
	if len(df):
		pvals = df["p"].to_numpy()
		if p_adjust == "fdr_bh":
			df["p_adj"] = _bh_fdr(pvals)
		elif p_adjust == "holm":
			df["p_adj"] = _holm(pvals)
		elif p_adjust == "none":
			df["p_adj"] = pvals
		else:
			raise ValueError("p_adjust must be 'none', 'holm', or 'fdr_bh'.")

		sig = (df["p_adj"] < 0.05).to_numpy()
		dirn = np.where(sig & (df["median_diff"] > 0), "↑ increase",
			np.where(sig & (df["median_diff"] < 0), "↓ decrease", "— no clear change"))
		df["direction"] = dirn

	# NEW: stash data & summaries for downstream plotting
	out = {
		"table": df,
		"baseline_index": int(baseline_index),
		"baseline_x": float(x[baseline_index]),
		"x_all": np.array(x, copy=True),
		"A_shape": A.shape,
		"A_row_mean": np.nanmean(A, axis=1),
		"A_row_median": np.nanmedian(A, axis=1),
		"n_mice": int(A.shape[1]),
		"k_levels": int(A.shape[0]),
	}
	if store_matrix:
		out["A"] = np.array(A, copy=True)  # optional, only if you really want the whole matrix back

	return out

# --- drop-in replacement: annotate_compare_to_baseline ---
import numpy as np
import matplotlib.pyplot as plt

def _stars_from_p(p):
    if not np.isfinite(p): return ""
    if p <= 0.001: return "***"
    if p <= 0.01:  return "**"
    if p <= 0.05:  return "*"
    return ""

def _format_p(p, p_fmt=None):
    """
    Format a p-value:
      - If p_fmt is a str like '.3g' or '.2e', uses that.
      - If None, uses auto: scientific for p<1e-3, else 3 significant figs.
    """
    if not np.isfinite(p):
        return "nan"
    if isinstance(p_fmt, str):
        return format(p, p_fmt)
    # auto
    return format(p, ".1e") if p < 1e-3 else format(p, ".3g")

def annotate_compare_to_baseline(
    ax,
    results_or_df,
    x=None,                 # optional: x for ALL levels; if None, uses results['x_all'] or df['x']
    y=None,                 # optional: y for ALL levels; if None, uses results['A_row_mean'] (fallback median)
    alpha=0.05,
    p_col="p_adj",
    show_ns=False,
    include_arrow_for_ns=False,
    offset_pts=(6, 4),
    fontsize=10,
    color="k",
    bbox=None,
    *,
    show_p=False,           # NEW: include p-value text in the label
    p_fmt=None              # NEW: formatting for p (e.g., '.3g', '.2e'); None => auto
):
    """
    Annotate asterisks and direction arrows for each non-baseline comparison.
    Now supports printing p-values (or q-values when p_col='p_adj') alongside the stars.

    Label examples:
      show_p=False           -> '** ↑'
      show_p=True, p_col='p' -> '** ↑ (p=0.004)'
      show_p=True, p_col='p_adj' -> '** ↑ (q=0.004)'
    """
    # Accept dict or DataFrame
    if isinstance(results_or_df, dict):
        results = results_or_df
        if "table" not in results:
            raise ValueError("Expected dict from compare_to_baseline(...), missing 'table'.")
        df = results["table"].copy()
    else:
        results = {}
        df = results_or_df.copy()

    if df.empty:
        return []

    xs_df   = df["x"].to_numpy()      if "x"   in df.columns else None
    rows_ix = df["row"].to_numpy()     if "row" in df.columns else None
    pvals   = df[p_col].to_numpy()     if p_col in df.columns else df["p"].to_numpy()
    meddiff = df["median_diff"].to_numpy()

    # Resolve full x and y if not provided
    if x is None:
        if "x_all" in results and rows_ix is not None:
            x_full = np.asarray(results["x_all"])
            xpos   = x_full[rows_ix]
        elif xs_df is not None:
            xpos = xs_df
        else:
            raise ValueError("Cannot infer x positions. Pass x=... or provide results with 'x_all'.")
    else:
        x = np.asarray(x)
        if rows_ix is not None and x.shape[0] > rows_ix.max():
            xpos = x[rows_ix]
        else:
            xpos = np.array([x[np.argmin(np.abs(x - xv))] for xv in xs_df])

    if y is None:
        if "A_row_mean" in results and rows_ix is not None:
            y_full = np.asarray(results["A_row_mean"], float)
            ypos   = y_full[rows_ix]
            va = "bottom"; y_offset_pts = abs(offset_pts[1])  # upward
        elif "A_row_median" in results and rows_ix is not None:
            y_full = np.asarray(results["A_row_median"], float)
            ypos   = y_full[rows_ix]
            va = "bottom"; y_offset_pts = abs(offset_pts[1])
        else:
            y_top = ax.get_ylim()[1]
            ypos  = np.full_like(xpos, y_top, dtype=float)
            va = "top"; y_offset_pts = -abs(offset_pts[1])    # downward
    else:
        y = np.asarray(y, float)
        if rows_ix is not None and y.shape[0] > rows_ix.max():
            ypos = y[rows_ix]
        else:
            if xs_df is None:
                raise ValueError("Cannot infer mapping to y. Pass y aligned to rows or provide df['x'].")
            ypos = np.array([y[np.argmin(np.abs(x - xv))] for xv in xs_df])
        va = "bottom"; y_offset_pts = abs(offset_pts[1])

    # Compose and draw labels
    artists = []
    p_name = "q" if p_col.lower() in ("p_adj", "q", "qval", "qvalue") else "p"
    for j in range(len(df)):
        p = pvals[j]
        stars = _stars_from_p(p)
        sig = (np.isfinite(p) and p <= alpha)

        # direction arrow
        if np.isfinite(meddiff[j]) and meddiff[j] > 0:
            arrow = "↑"
        elif np.isfinite(meddiff[j]) and meddiff[j] < 0:
            arrow = "↓"
        else:
            arrow = "•"

        # build label text
        if sig:
            label = f"{stars} {arrow}".strip()
        else:
            if show_ns:
                label = f"ns {arrow}".strip() if include_arrow_for_ns else "ns"
            else:
                if include_arrow_for_ns:
                    label = arrow
                else:
                    # skip entirely when not significant and user doesn't want arrows or ns
                    if not show_p:
                        continue
                    # If only showing p, we still need a base for the label; use arrow if available
                    label = arrow if arrow != "•" else ""

        # append p/q if requested
        if show_p and np.isfinite(p):
            label = f"{label} ({p_name}={_format_p(p, p_fmt)})".strip()

        txt = ax.annotate(
            label,
            xy=(xpos[j], ypos[j]),
            xytext=(offset_pts[0], y_offset_pts),
            textcoords="offset points",
            ha="left", va=va,
            fontsize=fontsize, color=color,
            bbox=bbox
        )
        artists.append(txt)

    return artists


# -------- Example usage --------
# results = utf.compare_to_baseline(A, x=w_zvert_list2use,
#                                   baseline_index=0,
#                                   test="wilcoxon",
#                                   p_adjust="holm",
#                                   store_matrix=False)  # you can set True if you really want A back
# utf.annotate_compare_to_baseline(
#     ax, results,
#     p_col="p_adj",
#     show_p=True,          # <—
#     p_fmt=".2e",          # show in scientific notation; try ".3g" for compact
#     show_ns=False,
#     include_arrow_for_ns=False,
#     offset_pts=(6, 6),
#     bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7)
# )

import numpy as np

def add_significance_stars(ax, x, y, p,
						alpha_levels=(0.05, 0.01, 0.001, 0.0001),
						symbols=('*', '**', '***', '****'),
						display='stars',
						add_legend=False,
						legend_loc='upper right',
						fontsize=14,
						x_offset=0.0,
						y_offset=0.02,
						color='black',
						ha='center',
						va='bottom'):

	import numpy as np

	x = np.asarray(x)
	y = np.asarray(y)
	p = np.asarray(p)

	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	dx = x_offset * (xlim[1] - xlim[0])
	dy = y_offset * (ylim[1] - ylim[0])

	for xi, yi, pi in zip(x, y, p):

		if np.isnan(pi):
			continue

		star = None
		for thresh, sym in reversed(list(zip(alpha_levels, symbols))):
			if pi < thresh:
				star = sym
				break

		text = None

		if display == 'stars':
			text = star

		elif display == 'pvalue':
			if pi < 0.001:
				text = "p<0.001"
			elif pi < 0.01:
				text = "p<0.01"
			elif pi < 0.05:
				text = "p<0.05"

		elif display == 'both':
			if star is not None:
				text = f"{star}\n(p={pi:.3g})"
			else:
				text = f"p={pi:.3g}"

		if text is not None:
			ax.text(xi + dx,
					yi + dy,
					text,
					fontsize=fontsize,
					color=color,
					ha=ha,
					va=va)

		# ---- Add small significance legend box (optional) ----
	if add_legend and display in ['stars', 'both']:

		legend_text = []
		for thresh, sym in zip(alpha_levels, symbols):
			legend_text.append(f"{sym}  p < {thresh}")

		legend_str = "\n".join(legend_text)

		loc_map = {
			'upper right': (0.98, 0.98, 'right', 'top'),
			'upper left':  (0.02, 0.98, 'left',  'top'),
			'lower right': (0.98, 0.02, 'right', 'bottom'),
			'lower left':  (0.02, 0.02, 'left',  'bottom'),
			'right':       (0.98, 0.5,  'right', 'center'),
			'left':        (0.02, 0.5,  'left',  'center'),
			'top':         (0.5,  0.98, 'center','top'),
			'bottom':      (0.5,  0.02, 'center','bottom')
		}

		if isinstance(legend_loc, str):
			xpos, ypos, ha_loc, va_loc = loc_map.get(
				legend_loc.lower(),
				(0.98, 0.98, 'right', 'top')
			)
		else:
			# custom tuple (x, y)
			xpos, ypos = legend_loc
			ha_loc, va_loc = 'left', 'bottom'

		ax.text(
			xpos, ypos,
			legend_str,
			transform=ax.transAxes,
			fontsize=fontsize * 0.8,
			ha=ha_loc,
			va=va_loc,
			bbox=dict(
				boxstyle="round,pad=0.4",
				facecolor="white",
				edgecolor="black",
				alpha=0.9
			),
			zorder=10
		)

##########################################################################