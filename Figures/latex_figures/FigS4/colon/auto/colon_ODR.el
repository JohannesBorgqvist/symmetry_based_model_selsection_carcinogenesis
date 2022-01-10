(TeX-add-style-hook
 "colon_ODR"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "crop" "convert={density=200,outext=.png}" "border=0.4cm")))
   (TeX-run-style-hooks
    "latex2e"
    "./Input/data"
    "./Input/IM-II_ODR"
    "./Input/PLM_ODR"
    "standalone"
    "standalone10"
    "pgfplots"
    "amsmath"
    "physics"
    "xcolor")
   (LaTeX-add-xcolor-definecolors
    "mixed_1"
    "mixed_2"
    "mixed_3"
    "pow_1"
    "pow_2"
    "pow_3"))
 :latex)

