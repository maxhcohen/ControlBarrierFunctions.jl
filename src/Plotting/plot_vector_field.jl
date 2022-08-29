"""
    plot_vector_field(
        xs,
        ys,
        f::Function;
        normalize_arrows=true,
        scale=0.15,
        lw=1,
        colormap=:viridis,
        colorbar=true,
    )

Plot the vector field `f(x,y)` over the region specified by xs and ys.
"""
function plot_vector_field(
    xs,
    ys,
    f::Function;
    normalize_arrows=true,
    scale=0.15,
    lw=1,
    colormap=:viridis,
    colorbar=true,
)
    # Make meshgrid from (xs, ys)
    Xs, Ys = meshgrid(xs, ys)

    # Compute vector field on the mesh
    dx, dy = mesh_vector_field(Xs, Ys, f, scale, normalize_arrows)

    # Map the magnitude of the vectorfield to a color
    c = vector_field_colors(Xs, Ys, f)
    
    return quiver(Xs, Ys, quiver=(dx, dy), line_z=c, lw=lw, c=colormap, colorbar=colorbar, arrows=true)
end

"""
    plot_vector_field!(
        xs,
        ys,
        f::Function;
        normalize_arrows=true,
        scale=0.15,
        lw=1,
        colormap=:viridis,
        colorbar=true,
    )

Plot the vector field `f(x,y)` over the region specified by xs and ys.
"""
function plot_vector_field!(
    xs,
    ys,
    f::Function;
    normalize_arrows=true,
    scale=0.15,
    lw=1,
    colormap=:viridis,
    colorbar=true,
)
    # Make meshgrid from (xs, ys)
    Xs, Ys = meshgrid(xs, ys)

    # Compute vector field on the mesh
    dx, dy = mesh_vector_field(Xs, Ys, f, scale, normalize_arrows)

    # Map the magnitude of the vectorfield to a color
    c = vector_field_colors(Xs, Ys, f)

    # Get the current figure
    fig = Plots.current()
    
    return quiver!(fig, Xs, Ys, quiver=(dx, dy), line_z=c, lw=lw, c=colormap, colorbar=colorbar)
end

"""
    mesh_vector_field(Xs, Ys, f::Function, scale::Float64, normalize_arrows::Bool)

Compute the vector field over the mesh specified by `Xs` and `Ys`.

If the norm of the vector field is zero anywhere then its likely we have some NaNs in dx and
dy. If so, just replace any `NaN` values with zero.
"""
function mesh_vector_field(Xs, Ys, f::Function, scale::Float64, normalize_arrows::Bool)
    # Check if we should normalize all vectors to unit length and then scale them up/down
    if normalize_arrows
        dfs = scale * normalize.(f.(Xs, Ys))
    else
        dfs = scale * f.(Xs, Ys)
    end

    # Reshape into individual vectors for the x and y components of vector field
    dx = [df[1] for df in dfs]
    dy = [df[2] for df in dfs]

    # Check if we divided by zero anywhere
    replace!(dx, NaN => 0.0)
    replace!(dy, NaN => 0.0)

    return dx, dy
end