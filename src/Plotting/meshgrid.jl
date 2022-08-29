"""
    meshgrid(xs, ys)

Create a meshgrid from the x and y coordinates specified by `xs`` and `ys`.`
"""
function meshgrid(xs, ys)
    Xs = [x for x in xs for y in ys]
    Ys = [y for x in xs for y in ys]

    return Xs, Ys
end