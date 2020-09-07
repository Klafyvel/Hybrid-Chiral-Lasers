% savefig.m
% A small utilitary to keep track of every figure I produce while being able to 
% fetch the latest ones quickly.
% 
% Save one fig in ./plots/[plot_type]/date/[plot_name].[ext] and one in ./plots/[plot_type]/latest/[plot_name].ext
% You can specify a format for saveas as last argument. By default it is `ext`.

function savefig(fig, plot_type, plot_name, ext, format)
    if nargin < 5
        format = ext;
    end
    warning('off', 'MATLAB:MKDIR:DirectoryExists')
    date_format = string(datetime, 'yyyy_MM_dd_hh_mm');
    mkdir('plots', plot_type);
    mkdir(sprintf('plots/%s', plot_type), date_format);
    mkdir(sprintf('plots/%s', plot_type), 'latest');
    warning('on', 'MATLAB:MKDIR:DirectoryExists')
    filename_1 = sprintf('plots/%s/%s/%s.%s', plot_type, date_format, plot_name, ext);
    filename_2 = sprintf('plots/%s/latest/%s.%s', plot_type, plot_name, ext);
    saveas(fig, filename_1, format)
    saveas(fig, filename_2, format)
end
