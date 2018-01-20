% Set general function for later application
function [fgen]=fgen(x1,x2)

fgen=exp(-pi.*(x1.^2+(x2-8).^2))+exp(-pi.*(x1.^2+(x2+8).^2));