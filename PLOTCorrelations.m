function [rho,pval] = PLOTCorrelations(X,Y,CorrType,Color)

s = scatter(X,Y,10,'filled');
s.LineWidth= 1;
if nargin == 4
        s.MarkerFaceColor = Color;
else
    s.MarkerFaceColor = 'k';
end

if nargin >= 3
    [rho,pval] = corr(X,Y,'Type',CorrType);
else
    [rho,pval] = corr(X,Y);
end
rho = round(rho,2);
hold on
h = lsline;
h.LineWidth = 2;
if pval < 0.05
    SigText = strcat(num2str(rho),'^{*}');
    if pval < 0.01
        SigText = strcat(num2str(rho),'^{**}');
        if pval < 0.01
            SigText = strcat(num2str(rho),'^{***}');
        end
    end
else
    SigText = 'n.s.';
end
CoordinatesY = ylim;
CoordinatesY = CoordinatesY(1,1) + (CoordinatesY(1,2)-CoordinatesY(1,1))*0.1;
CoordinatesX = xlim;
CoordinatesX = CoordinatesX(1,1) + (CoordinatesX(1,2)-CoordinatesX(1,1))*0.5;
text(CoordinatesX,CoordinatesY,strcat(SigText,num2str(round(pval,2))));