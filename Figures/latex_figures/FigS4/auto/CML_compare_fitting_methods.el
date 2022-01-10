(TeX-add-style-hook
 "CML_compare_fitting_methods"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("standalone" "tikz" "crop" "convert={density=300,outext=.png}" "border=0.3cm" "width=18cm" "height=3cm")))
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

