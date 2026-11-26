"""
Build the documentation as a pdf using LaTeX

Note this requires LuaTeX to be installed, e.g. on Ubuntu or Mint
```
sudo apt install texlive-luatex
```
It may also need the `texlive-latex-extra` [JTO: already had this installed, so haven't
tested without it].
"""

using Pkg

repo_dir = dirname(dirname(@__FILE__))
Pkg.develop([PackageSpec(path=joinpath(repo_dir, "moment_kinetics")),
             PackageSpec(path=joinpath(repo_dir, "makie_post_processing", "makie_post_processing")),
             PackageSpec(path=joinpath(repo_dir, "plots_post_processing", "plots_post_processing"))])
Pkg.instantiate()

using Documenter
using moment_kinetics, makie_post_processing, plots_post_processing
include("../util/MKJacobianUtils.jl")
include("organise_pages.jl")

if get(ENV, "CI", nothing) == "true"
    latex_kwargs = (platform = "docker",)
else
    latex_kwargs = ()
end

makedocs(
    sitename = "KfT",
    format = Documenter.LaTeX(; latex_kwargs...),
    modules = [moment_kinetics, makie_post_processing, plots_post_processing, MKJacobianUtils],
    authors = "KfT contributors",
    pages = organise_pages()
)
