@doc raw"""
    SmithMethod

Warren D. Smith (Jan 1989) enumeration method.

Simplified Steiner tree backtrack code.

SMT representation: Steiner point i is adjacent to points adj[i][0..2].
 * (If neighbor is Steiner - p added!)
 * There are N-2 Steiner points and N regular sites.
 * The coordinates of Steiner point i are XX[i+p][0..N-1], and
 * The coordinates of regular site i are XX[i][N-1], i=1,2,..
 * Edge i has endpoints edge[i][0] < edge[i][1] (sam p-add convention).
""" struct SmithMethod <: ExactMethod end

const ε₀ = 0.0001
const ε₁ = 0.005
const ε₂ = 0.001
const ε₃ = 0

mutable struct SmithData{T,U <: Integer}
	# Dimensions
	k::Int # ?
	N::Int # space dimension
	p::Int # terminal points
	s::Int # Steiner points
	n::Int # total points

	# Coordinates
	x::Matrix{T} # points

	# Topology
	θ::Vector{U} # topology
	ψ::Vector{U} # best vector
	ℵ::Matrix{U} # adjacency
	δ::Vector{U} # valence
	E::Matrix{U} # edges
	S::Vector{U} # stack

	# Ancillary
	u::Vector{U}
	v::Matrix{T}
	B::Matrix{T}
	C::Matrix{T}

	# Lengths
	L::Matrix{T}
	Δ::Vector{T}
	ℓ::T # any upper bound on the length of the SMT

	# Constants
	α::T # Scale
	β::T # ?

	# Random
	rng::Random.AbstractRNG

	function SmithData{T,U}(points::Matrix{T}; rng::Random.AbstractRNG = Random.Xoshiro()) where {T,U <: Integer}
		# Dimensions
		k = 3
		N = size(points, 1) # space dimension
		p = size(points, 2) # terminal points
		s = p - 2           # Steiner points
		n = p + s           # total points

		# Coordinates
		x = Matrix{T}(undef, N, 2n)

		# Store coordinates for terminal points
		x[1:N, 1:p] .= points[1:N, 1:p]

		# Ancillary
		u = Vector{U}(undef, 3)
		v = Matrix{T}(undef, N, 3)
		B = Matrix{T}(undef, 3, n)
		C = Matrix{T}(undef, N, n)

		# Topology
		ψ = Vector{U}(undef, n)        # best topology vector
		θ = Vector{U}(undef, n)        # topology vector
		ℵ = Matrix{U}(undef, 3, n - 2) # adjacency list
		δ = Vector{U}(undef, s)        # Valence
		E = Matrix{U}(undef, 2, 2n)    # edges
		S = Vector{U}(undef, n ^ 2)    # stack
		
		# Lengths
		L = Matrix{T}(undef, 3, n)     # edge length
		Δ = Vector{T}(undef, 3)        # ?
		ℓ = typemax(T)                 # length upper bound

		# Constants
		α = zero(T) # Scale
		β = zero(T) # ?

		for i = 1:p, j = 1:N
			α = max(α, abs(x[j, i] - x[j, 1]))
		end

		return new{T,U}(
			# Dimensions
			k, N, p, s, n,

			# Coordinates
			x,

			# Topology
			θ, ψ, ℵ, δ, E, S,

			# Ancillary
			u, v, B, C,


			# Lengths
			L, Δ, ℓ,

			# Constants
			α, β,

			# Random
			rng,
		)
	end
end

SmithData(args...; kw...) = SmithData{Float64,Int}(args...; kw...)

function smith_solve_3!(s::SmithData{T,U}) where {T,U}
	smith_build_tree!(s, 0)

	@show ℓ = smith_length!(s)
	@show ε = smith_error(s)

	while true
		@show τ = ε₀ * ε / s.p

		smith_optimize!(s, τ)

		@show ℓ = smith_length!(s)
		@show ε = smith_error(s)

		ε > ℓ * ε₀ || break
	end

	return nothing
end

function smith_solve!(s::SmithData{T,U}) where {T,U}
	# special case
	s.p == 3 && return smith_solve_3!(s)

	error("p ≤ 3, please")
	
	# HINT: optionally, sort sites in some nice order here

	# # constants
	# k = 1
	# m = 0
	
	# while true
	# 	nc = 0

	# 	l = 2k + 1

	# 	while l > 1 && k <= s.p - 3
	# 		s.θ[k] = l

	# 		# Build the tree represented by the topol. vector θ[1..k]
	# 		smith_build_tree!(s, k)

	# 		@label(iter)

	# 		# .. and optimize it until either obviously bad or small error figure happens
	# 		ℓ = smith_length!(s)
	# 		ε = smith_error(s)

	# 		if (ℓ - s.ℓ) < ε
	# 			if ε₁ * ℓ < ε
	# 				smith_optimize!(s, ε₀ * ε / s.p)

	# 				@goto(iter)
	# 			end

	# 			if k >= s.p - 3
	# 				while true
	# 					smith_optimize!(s, ε₀ * ε / s.p)

	# 					ℓ = smith_length!(s)
	# 					ε = smith_error(s)

	# 					ε > ℓ * ε₀ || break
	# 				end

	# 				if ℓ < s.ℓ
	# 					# update best solution
	# 					s.ℓ = ℓ 
	# 					s.ψ .= s.θ
	# 				else
	# 					i = nc
	# 					s.k += 1

	# 					while i > 1 && s.L[i] < ℓ
	# 						s.S[m + i + 1] = s.S[m + i]
	# 						s.L[i + 1]     = s.L[i]

	# 						i -= 1
	# 					end

	# 					i += 1
	# 					s.S[m + i] = l
	# 					s.L[i] = ℓ
	# 				end
	# 			end
	# 		end
	# 		l -= 1
	# 	end

	# 	m += nc

	# 	while nc <= 0
	# 		k -= 1
	# 		k <= 0 && return nothing

	# 		nc = s.S[m]
	# 		m -= 1
	# 	end

	# 	s.A[k] = s.S[m]
	# 	s.S[m] = nc - 1

	# 	if k < p - 3
	# 		k += 1
	# 	else
	# 		m -= 1
	# 	end
	# end

	return nothing
end

function smith_prep!(s::SmithData, i::Integer, k::Integer)
	l = s.u[k]
	δ = s.Δ[k]

	if (l > s.p)
		# if l points to a Steiner point
		s.δ[i]  += 1 # increase the valence for `i`
		s.B[k,i] = δ #
	else
		# l points to a terminal point
		s.C[:,i] .+= s.x[:,l] .* δ
	end
end

function smith_optimize!(s::SmithData{T,U}, τ::Float64) where {T,U}
	# τ: a small positive number */
	# finds better coordinates XX[p+1..p+k1][] for the k1 Steiner points
	# of tree T by: doing a relaxation iteration. Assumes that edge lengths of old tree
	# have been per-stored in array s.L[][] */

	K = trunc(U, s.β - 2)

	stack = sizehint!(U[], s.n)
	queue = sizehint!(U[], s.n)

	return nothing

	# First: compute B array, C array, and valences. Set up leafQ.
	for i = K:-1:1
		s.u[:] = s.ℵ[:,i] # adjacency list for `i`
		s.Δ[:] = 1 ./ (s.L[:,i] .+ τ)

		# Have now figured out reciprocal distances q0, q1, q2 from
		# Steiner pt. i to neighbors n0, n1, n2 */
		s.Δ[:] ./= sum(s.Δ)

		s.δ[i] = zero(U)

		s.B[:,i] .= zero(T)
		s.C[:,i] .= zero(T)

		smith_prep!(s, i, 1)
		smith_prep!(s, i, 2)
		smith_prep!(s, i, 3)

		# Now: Steiner point i has Steiner valence val[i];
		# coords obey eqns
		# XX[i+p][] = sum(j) of B[j,i] * XX[nj][] + C[i][] */
		if (s.δ[i] <= 1)
			push!(queue, i)
		end
	end

	# Have set up equations - now to solve them.
	# Second: eliminate leaves
	while !isempty(queue)
		# pick leaf vertex
		iₜ = popfirst!(queue)

		s.δ[iₜ] -= 1 # reduce valence for `iₜ`

		# get related Steiner point
		iₛ = iₜ + s.p

		# now to eliminate leaf `iₜ`
		push!(stack, iₜ) # push `iₜ` onto stack
		
		j = something(findfirst(!iszero, s.B[:, iₜ]), 3)

		s.Δ[1] = s.B[j,iₜ]
		
		j = 1 + (s.ℵ[j,iₜ] - s.p) # neighbor is j

		s.δ[j] -= 1 # reduce valence for `j`

		if s.δ[j] == 1
			push!(queue, j)
		end # new leaf?

		m = something(findfirst(==(iₛ), s.ℵ[:,j]), 3)

		s.Δ[2] = s.B[m,j]

		s.B[m,j] = 0.0

		γ = 1.0 / (1.0 - s.Δ[2] * s.Δ[1])		
		
		s.B[:,j] .= γ .* (s.B[:,j])
		s.C[:,j] .= γ .* (s.C[:,j] .+ s.Δ[2] .* s.C[:,iₜ])
	end

	# III: Solve 1-vertex tree!
	iₜ = queue[begin]
	iₛ = iₜ + s.p

	s.x[:,iₛ] .= s.C[:,iₜ]

	# IV: backsolve
	while !isempty(stack)
		iₜ = pop!(stack)
		iₛ = iₜ + s.p

		j = something(findfirst(!iszero, s.B[:, iₜ]), 3)

		s.Δ[1] = s.B[j,iₜ]

		j = s.ℵ[j,iₜ] # neighbor is j */

		s.x[:,iₛ] .= s.C[:,iₜ] + s.Δ[1] * s.x[:,j]
	end

	return
end

function smith_build_tree!(s::SmithData{T,U}, k::Integer) where {T,U}
    # builds tree represented by topvec[1..k]. Initial location of new Steiner pts
    # is a slightly random perturbation of the centroid of its neighbors
	
    # First build the tree corresponding to the null vector
	
	# ! Connect the first 3 terminal points to a Steiner point
	# ! chosen to be a perturbation away from the center
	s.β = 3
	m   = s.p + 1 #

	s.ℵ[1:3,1] .= 1:3
	s.E[1,1:3] .= 1:3
	s.E[2,1:3] .= m

	for i = 1:s.N
		s.x[i,m] = sum(s.x[i,1:3]) / 3.0 + ε₂ * s.α * rand(s.rng)
	end

	for i = 1:k
		# Now: do vector element topvec[i]
		en = i + 3
		m  = i + 1
		sn = m + s.p
		e  = s.θ[i]
		ea = s.E[1,e]
		eb = s.E[2,e]
		s.ℵ[1,m] = ea
		s.ℵ[2,m] = eb
		s.ℵ[3,m] = en
		m = ea - s.p

		if m > 0
			for j = 1:3
				if s.ℵ[j,m] == eb
					s.ℵ[j,m] = sn
					break
				end
			end
		end

		m = eb - s.p

		if m > 0
			for j = 1:3
				if s.ℵ[j,m] == ea
					s.ℵ[j,m] = sn
					break
				end
			end
		end

		s.E[2,e] = sn
		
		e = 2en - 4

		s.E[1,e] = en
		s.E[2,e] = sn

		e += 1

		s.E[1,e] = eb
		s.E[2,e] = sn

		for j = 1:s.N
			s.x[j,sn] = (s.x[j,ea] + s.x[j,eb] + s.x[j,en]) / 3.0 + ε₂ * s.α * rand(s.rng)
		end
	end

	s.β = k + 3
	
	# Tree is now built. Initial coords in general position.
	
	return
end

function smith_dist(s::SmithData{T,U}, i::Integer, j::Integer) where {T,U}
	return norm(s.x[:, i] - s.x[:, j])
end

function smith_length_delta!(s::SmithData{T,U}, m::Integer, i::Integer, k::Integer) where {T,U}
	l = s.u[m]

	if l < k
		Δℓ = smith_dist(s, l, k)
		
		s.L[m,i] = Δℓ

		l -= s.p

		if l > 0
			for j = 1:3
				if s.ℵ[j,l] == k
					s.L[j,l] = Δℓ
					break
				end
			end
		end

		s.u[m] = l

		return Δℓ
	else
		return zero(T)
	end
end

function smith_length!(s::SmithData{T,U}) where {T,U}
	# Stores edge lengths of tree T in array s.L[1..k1][0..2] and returns total length.
	m = trunc(U, s.β - 2)
	ℓ = zero(T)

	for i = 1:m
		s.u[:] .= s.ℵ[:, i]

		ℓ += smith_length_delta!(s, m, i, i + s.p)
	end
	
	# Have now figured out distance s.L[i][00.3] from Steiner pt. i to neighbors. */
	
	return ℓ
end

function smith_error(s::SmithData{T,U}) where {T,U}
	# Returns the error figure of tree T with Steiner coords in XX[][].
	# Assumes edge lengths have been pre-stored in array s.L[][]. */
	m = s.k - 2
	ε = zero(T)

	for i = 1:m
		k = i + s.p

		s.u[:]   .= s.ℵ[:,i]
		s.v[:,:] .= s.x[:,s.u] .- s.x[:,k]

		for (l₁, l₂) in [(1,2), (1,3), (2, 3)]
			let Δε = 2 * (s.v[:,l₁]' * s.v[:,l₂]) + s.L[i,l₁] * s.L[i,l₂]
				# only angles < 120 cause error */
				ε += max(zero(T), Δε)
			end
		end
	end
	
	return √ε
end

function solve(::SmithMethod, points::Matrix{K}) where {K}
	s = SmithData{K,Int}(points; rng = Random.Xoshiro())

	smith_solve!(s)

	return s # Return the steiner tree here! 
end
