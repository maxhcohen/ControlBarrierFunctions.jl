"""
	HighOrderCBF

Type representing a High Order Control Barrier Function (HOCBF).

This functionality allows for automatically computing higher order derivatives of a given constraint function `h` to extend it to an HOCBF.

# Fields
- `h::Function` : function defining the constraint set `h(x) ≥ 0`
- `relative_degree::Int` : relative degree of the constraint function
- `cbf::ControlBarrierFunction` : CBF used to represent HOCBF constraint (which is not actually a CBF)
"""
struct HighOrderCBF
	h::Function
	relative_degree::Int
	cbf::ControlBarrierFunction
end


###### HighOrderCBF constructors #####

"""
	HighOrderCBF(h::Function, Σ::ControlAffineSystem, relative_degree::Int, α::Vector{Float64})

Constructs a high order CBF with a given relative degree and collection of linear class K functions.
"""
function HighOrderCBF(
	h::Function, Σ::ControlAffineSystem, relative_degree::Int, α::Vector{Float64},
)

	# Make sure we have enough class K functions
	@assert relative_degree == length(α) "Didn't pass in enough class K functions!"

	# Allocate a vector of functions and initialize first component with h
	arr = Array{Function}(undef, relative_degree)
	arr[1] = h

	# Loop through array to recursively define HOCBF
	for i in 1:(relative_degree-1)
		h_old = arr[i]
		h_new = x -> ForwardDiff.gradient(h_old, x)' * Σ.f(x) + α[i] * h_old(x)
		arr[i+1] = h_new
	end

	# Make CBF from last function in array
	cbf = ControlBarrierFunction(arr[relative_degree], Σ, r -> α[relative_degree] * r)

	return HighOrderCBF(h, relative_degree, cbf)
end

"""
	HighOrderCBF(h::Function, Σ::ControlAffineSystem, relative_degree::Int, α::Function)

Constructs a high order CBF with a given relative degree and collection of general class K functions.

Note that the `α` function must be vector-valued, e.g., `α(r) = [α1(r), α2(r)]`, where `α1` and `α2` class K functions.  
"""
function HighOrderCBF(
	h::Function, Σ::ControlAffineSystem, relative_degree::Int, α::Function,
)
	# Allocate a vector of functions and initialize first component with h
	arr = Array{Function}(undef, relative_degree)
	arr[1] = h

	# Loop through array to recursively define HOCBF
	for i in 1:(relative_degree-1)
		h_old = arr[i]
		h_new = x -> ForwardDiff.gradient(h_old, x)' * Σ.f(x) + α(h_old(x))[i]
		arr[i+1] = h_new
	end

	# Make CBF from last function in array
	cbf = ControlBarrierFunction(arr[relative_degree], Σ, r -> α(r)[relative_degree])

	return HighOrderCBF(h, relative_degree, cbf)
end
