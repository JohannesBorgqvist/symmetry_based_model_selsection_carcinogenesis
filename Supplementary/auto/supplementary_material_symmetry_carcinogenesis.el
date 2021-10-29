(TeX-add-style-hook
 "supplementary_material_symmetry_carcinogenesis"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "12pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "margin=1in") ("parskip" "parfill") ("inputenc" "utf8") ("fontenc" "T1") ("helvet" "scaled=1") ("sfmath" "helvet") ("footmisc" "symbol")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "./Input/derivation_symmetries"
    "./Input/data_processing"
    "./Input/num_val"
    "article"
    "art12"
    "geometry"
    "parskip"
    "inputenc"
    "fontenc"
    "helvet"
    "amsmath"
    "amssymb"
    "amsfonts"
    "amsthm"
    "mathtools"
    "sfmath"
    "physics"
    "footmisc"
    "url"
    "graphicx"
    "sansmathfonts"
    "sansmath")
   (LaTeX-add-labels
    "Oxford"
    "WennerGren"
    "Linacre")
   (LaTeX-add-bibliographies
    "model_selection_symmetries_cancer"))
 :latex)

