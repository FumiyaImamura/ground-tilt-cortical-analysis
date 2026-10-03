classdef plotCartInPolar <matlab.mixin.Copyable
	properties
		ax_polar
		ax_cart
		hlinks
	end
	methods
		function obj=plotCartInPolar()
			obj.hlinks=cell(0,1);
		end
		function makeCartAxes(obj,ax_polar)
			arguments
				obj
				ax_polar=[]
			end
			if isempty(ax_polar)
				ax_polar=gca();
				if ~isequal(ax_polar.Type,'polaraxes')
					error("The current axes need to be polar")
				end
			end
			obj.ax_polar=ax_polar;
			ax_polar=obj.ax_polar;
			ax_cart = axes();
			ax_cart.Position = ax_polar.Position;
			axes(ax_cart);
			obj.ax_cart=ax_cart;
			obj.hlinks{end+1}=linkprop([ax_polar,ax_cart],"Position");
		end
		function [x,y]=pol2cart(obj,theta,rho)
			arguments
				obj
				theta
				rho
			end
			rho=rho-obj.ax_polar.RLim(1);
			[x,y]=pol2cart(theta,rho);
		end
		function closeCartAxes(obj)
			arguments
				obj
			end
			ax_cart=obj.ax_cart;
			ax_polar=obj.ax_polar;

			axes(ax_cart);
			ax_cart.Visible="off";
			Rspan=ax_polar.RLim(2)-ax_polar.RLim(1);
			ylim([-Rspan,Rspan])
			xlim([-Rspan,Rspan])
			axis square
			axes(ax_polar)
			uistack(obj.ax_cart,'top')
		end
	end
end

function UnitTest()
%% Unit test
figure;polarplot([0,1,2,3,4],[0,1,2,3,4],LineWidth=2)
rlim([0.5,4])
a=plotCartInPolar();
a.makeCartAxes()
[x,y]=a.pol2cart([0,1,2,3,4],[0,1,2,3,4]);
scatter(x,y,"filled")
a.closeCartAxes()
%% test on tiledlayout
figure(WindowState="maximized");tiledlayout('flow')
nexttile;polarplot([0,1,2,3,4],[0,1,2,3,4],LineWidth=2)
rlim([0,4])
a=plotCartInPolar();
a.makeCartAxes()
[x,y]=a.pol2cart([0,1,2,3,4],[0,1,2,3,4]);
scatter(x,y,"filled");hold on
plot([0,1],[0,2])
a.closeCartAxes()

nexttile();
theta=[1:32]/8;
rho=2+rand(32,8);
polarplot(theta,rho)
a.makeCartAxes(gca);
[Mean,SE]=calculateMeanSE(rho,1,x=theta);
[mean_x,mean_y]=a.pol2cart(Mean.x,Mean.y);
[se_x,se_y]=a.pol2cart(SE.x,SE.y);
plot(mean_x,mean_y,LineWidth=2)
hold on;
patch(se_x,se_y,"red",facealpha=0.3,edgealpha=0)
a.closeCartAxes()
%%


end
