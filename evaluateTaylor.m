TL1=evalin(symengine,'mtaylor(sin(x^2+y),[x,y],8)')
% I am trying to run the function mtaylor from the MuPAD engine in MatLab,
% which provides a multivariate Taylor expansion of a function.字符串类型计算结果.
Fxy=sym('sin(x^2 + y)')
Fxy_TL1 = Fxy - TL1
figure(1)
ezsurf(Fxy,[-2,2,-3,3])
shading interp
% The shading function controls the color shading of surface
% and patch graphics objects.
% shading interp varies the color in each line segment and face by
% interpolating the colormap index or true color value across the line or
% face.
view([-63,52])
colormap(spring)
% colormap name sets the colormap for the current figure to the built-in
% colormap specified by name.
light,light('position',[-10,4,50],'style','local','color','r')
% light('PropertyName',propertyvalue,...)
% light creates a light object in the current axes. Lights affect only
% patch and surface objects.
% 
% light('PropertyName',propertyvalue,...) creates a light object using the
% specified values for the named properties. For a description of the
% properties, see Light Properties. The MATLAB? software parents the light
% to the current axes unless you specify another axes with the Parent
% property.
% 
% handle = light(...) returns the handle of the light object created.
figure(2)
ezsurf(TL1,[-2,2,-3,3])
shading interp
view([-43,54])
colormap(spring)
light
light('position',[-10,2,2],'style','local','color',[0.8,0.3,0.3])
light('position',[-2,-10,2],'style','local','color',[0.4,0.5,0.7])

figure(3)
ezsurf(Fxy,[-0.5,0.5,-0.5,0.5],'circ')
axis([-1,1,-1,1,-2,2])
shading interp
colormap(spring)
view([-49,17])
light
light('position',[-30,0,-2],'style','local','color','r')

figure(4)
ezsurf(Fxy_TL1,[-0.5,0.5],'circ')
shading interp
colormap(spring)
view([-53,34])
light
light('position',[-10,15,0],'style','local','color',[0.8,0.3,0.3])
% In calculus, Taylor's theorem gives an approximation of a k-times
% differentiable function around a given point by a k-th order Taylor
% polynomial. For analytic functions the Taylor polynomials at a given
% point are finite order truncations of its Taylor series, which completely
% determines the function in some neighborhood of the point. The exact
% content of "Taylor's theorem" is not universally agreed upon. Indeed,
% there are several versions of it applicable in different situations, and
% some of them contain explicit estimates on the approximation error of the
% function by its Taylor polynomial. Taylor's theorem is named after the
% mathematician Brook Taylor, who stated a version of it in 1712. Yet an
% explicit expression of the error was not provided until much later on by
% Joseph-Louis Lagrange. An earlier version of the result was already
% mentioned in 1671 by James Gregory.[1] Taylor's theorem is taught in
% introductory level calculus courses and it is one of the central
% elementary tools in mathematical analysis. Within pure mathematics it is
% the starting point of more advanced asymptotic analysis, and it is
% commonly used in more applied fields of numerics as well as in
% mathematical physics. Taylor's theorem also generalizes to multivariate
% and vector valued functions f\colon \mathbb R^n\rightarrow\mathbb R^m on
% any dimensions n and m. This generalization of Taylor's theorem is the
% basis for the definition of so-called jets which appear in differential
% geometry and partial differential equations.