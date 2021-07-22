(TeX-add-style-hook
 "supplementary_material_symmetry_carcinogenesis"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=0.7in") ("parskip" "parfill") ("inputenc" "utf8")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "geometry"
    "parskip"
    "inputenc"
    "amsmath"
    "amssymb"
    "amsfonts"
    "amsthm"))
 :latex)

