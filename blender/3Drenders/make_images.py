import subprocess
from slugify import slugify

LATEX = "pdflatex"

TEMPLATE = r"""
\documentclass[border=1pt]{standalone}
\usepackage{lmodern} 
\usepackage[intlimits]{amsmath}
\usepackage{amsthm, amssymb, amsfonts} 
\begin{document}
%s
\end{document}
"""

expressions = (
        "$n_1$",
        "$n_2$",
        r"$\bar{n},\;\delta n$",
        r"$\bar{n},\;\delta n,\;\psi$",
        r"$L = L_p\times m$",
        "$L_p/2$",
        "$n_L$",
        "$n_R$",
        r"$\bar{n}_L,\;\delta n_L$",
        r"$\bar{n}_R,\;\delta n_R,\;\psi$",
        "$L_L$",
        "$L_R$",
        "defect",
        "$\mathbf{\hat{x}}$",
        "$\mathbf{\hat{y}}$",
        "$\mathbf{\hat{z}}$",
        "$n_o$",
        "$n_e$"
)

for expr in expressions:
    print("Doing", expr)
    latex = subprocess.Popen([LATEX], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    latex.communicate((TEMPLATE % expr).encode('utf8'))
    cairo =  subprocess.Popen(["pdftocairo", "-svg", "texput.pdf", slugify(expr, separator='_')+'.svg'])

