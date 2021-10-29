(TeX-add-style-hook
 "FigS6"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "border=0.3cm")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "left=2.2cm" "right=2.2cm" "top=2.5cm" "bottom=2.0cm" "a4paper")))
   (TeX-run-style-hooks
    "latex2e"
    "standalone"
    "standalone10"
    "geometry"
    "pgfplots"
    "amsmath"
    "physics"))
 :latex)

