using Glob

"""
    organise_pages()

Organise the files in `src/` into a nested list so that the sidebar of the docs page has a
nice structure.
"""
function organise_pages()
    doc_files = glob("src/*.md")

    # Remove the src/ prefix
    doc_files = [basename(s) for s ∈ doc_files]

    function create_section(name, section_pages::Vector)
        # Check all of `section_pages` actually exist.
        if !all(p ∈ doc_files for p ∈ section_pages)
            error("Some of the pages $section_pages were not present in the existing doc "
                  * "files $doc_files.")
        end

        # Remove section_pages from doc_files, so that we can do something with all the
        # entries that are not explicitly handled.
        filter!(x -> x ∉ section_pages, doc_files)

        if length(section_pages) == 1
            return name => section_pages[1]
        else
            return name => section_pages
        end
    end

    pages = Any[]

    push!(pages, create_section("Home", ["index.md"]))

    push!(pages, create_section("Getting started", ["getting_started.md",
                                                    "machine_setup_notes.md",
                                                    "known_issues.md"]))

    push!(pages, create_section("Post processing", ["post_processing_notes.md"]))

    push!(pages, create_section("Equations", ["moment_kinetic_equations.md",
                                              "moment_kinetic_equations_electrons.md"]))

    developers = create_section("Developers", ["developer-index.md",
                                               "developer-intro.md",
                                               "debugging-hints.md",
                                               "manufactured_solution_test_examples.md",
                                              ])

    modules = create_section("Modules", filter(x -> startswith(x, "zz_"), doc_files))

    for p ∈ doc_files
        push!(pages, p)
    end

    push!(pages, developers)
    push!(pages, modules)

    return pages
end
