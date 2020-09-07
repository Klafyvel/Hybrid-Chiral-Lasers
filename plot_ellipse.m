% plot_ellipse.m
% 
% Plot a mode given theta, epsilon, polarization and its position.

function plot_ellipse(epsilon, theta, lcp, x, y,c,marker,hide_text)
    if nargin < 8
        hide_text = false;
    end
    for j=1:length(epsilon)
        if lcp(j)
            s = sprintf('$\\leftarrow$ %d (LCP)', j);
            if nargin < 6
                c = 'red';
            end
            ls = ':';
            m = '<';
        else
            s = sprintf('$\\leftarrow$ %d (RCP)', j);
            if nargin < 6
                c = 'blue';
                ls = ':';
            else
                ls='-';
            end
            m = '>';
        end
        [x_ellipse,y_ellipse] = ellipse(epsilon(j), theta(j), x(j), y(j));
        plot(x_ellipse, y_ellipse,'Color', c, 'LineStyle', ls,'LineWidth',2)
        hold on
        [max_y, max_x_index] = max(y_ellipse);
        max_x = x_ellipse(max_x_index);
        scatter(max_x,max_y,m,'filled','MarkerFaceColor',c)
        labels{j} = {s};
    end
    if nargin > 6
        scatter(x,y,marker,'filled','MarkerEdgeColor',c,'MarkerFaceColor',c)
    elseif nargin == 6
        scatter(x,y,'filled','MarkerEdgeColor',c,'MarkerFaceColor',c)
    else
        scatter(x,y,'filled')
    end
    if ~hide_text
        text(x,y,labels, 'interpreter', 'latex')
    end
end
