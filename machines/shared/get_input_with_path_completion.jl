"""
    get_input_with_path_completion(message=nothing)

Print `message` and get user input using the `read` utility from bash to allow
path-completion.

Solution adapted from
https://discourse.julialang.org/t/collecting-all-output-from-shell-commands/15592/2
"""
function get_input_with_path_completion(message=nothing)
    if message !== nothing
        println(message)
    end

    # The `out` object is used by `pipeline()` to capture output from the shell command.
    out = Pipe()

    # Use Julia's shell command functionality to actually run bash - not sure how to use
    # `read` and get output from it without using bash.
    run(pipeline(`bash -c "read -e -p '> ' USERINPUT; echo \$USERINPUT"`, stdout=out))

    # Need to close `out.in` to be able to read from `out`.
    close(out.in)

    # chomp removes the trailing '\n'. 'String()' converts Vector{Char} to a String.
    input = chomp(String(read(out)))

    return input
end
