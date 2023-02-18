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

struct SmithData{T,U <: Integer}
	# Dimensions
	k::Int
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
	E::Matrix{U} # edges
	S::Vector{U} # stack

	# Ancillary
	u::Vector{U}
	v::Matrix{T}
	B::Matrix{T}
	C::Matrix{T}

	# Constants
	α::T # Scale

	# Lengths
	L::Matrix{T}
	ℓ::T # any upper bound on the length of the SMT

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
		
		α = zero(T) # Scale

		for i = 1:p, j = 1:N
			α = max(α, abs(x[j, i] - x[j, 1]))
		end

		# Ancillary
		u = Vector{U}(undef, 3)
		v = Matrix{T}(undef, N, 3)
		B = Matrix{T}(undef, 3, n)
		C = Matrix{T}(undef, N, n)

		# Topology
		ψ = Vector{U}(undef, n)
		θ = Vector{U}(undef, n)
		ℵ = Matrix{U}(undef, 3, n - 2)
		E = Matrix{U}(undef, 2, 2n)
		S = Vector{U}(undef, n ^ 2)
		
		# Lengths
		L = Matrix{T}(undef, 3, n)
		ℓ = typemax(T)

		return new{T,U}(
			# Dimensions
			k, N, p, s, n,

			# Coordinates
			x,

			# Topology
			θ, ψ, ℵ, E, S,

			# Ancillary
			u, v, B, C,

			# Constants
			α,

			# Lengths
			L, ℓ,

			# Random
			rng,
		)
	end
end

SmithData(args...; kw...) = SmithData{Float64,Int}(args...; kw...)

function smith_solve_3!(s::SmithData{T,U}) where {T,U}
	smith_build_tree!(s, 0)

	ℓ = smith_length!(s)
	ε = smith_error(s)

	while true
		smith_optimize!(s, ε₀ * ε / s.p)

		ℓ = smith_length!(s)
		ε = smith_error(s)

		ε > ℓ * ε₀ || break
	end

	return nothing
end

function smith_solve!(s::SmithData{T,U}) where {T,U}
	# special case
	s.p == 3 && return smith_solve_3!(s)
	
	# HINT: optionally, sort sites in some nice order here

	# constants
	δ = 

	k = 1
	m = 0
	
	while true
		nc = 0

		l = 2k + 1

		while l > 0 && k <= s.p - 3
			s.θ[k] = l

			# Build the tree represented by the topol. vector θ[1..k]
			smith_build_tree!(s, k)

			@label(iter)

			# .. and optimize it until either obviously bad or small error figure happens
			ℓ = smith_length!(s)
			ε = smith_error(s)

			if (ℓ - s.ℓ) < ε
				if ε₁ * ℓ < ε
					smith_optimize!(s, ε₀ * ε / p)

					@goto(iter)
				end

				if k >= p - 3
					while true
						smith_optimize!(s, ε₀ * ε / p)

						ℓ = smith_length!(s)
						ε = smith_error(s)

						ε > ℓ * ε₀ || break
					end

					if ℓ < s.ℓ
						# update best solution
						s.ℓ = ℓ 
						s.ψ .= s.θ
					else
						i = nc
						n += 1

						while i > 1 && s.L[i] < ℓ
							s.S[m + i + 1] = s.S[m + i]
							s.L[i + 1]     = s.L[i]

							i -= 1
						end

						i += 1
						s.S[m + i] = l
						s.L[i] = ℓ
					end
				end
			end
			l -= 1
		end

		m += nc

		while nc <= 0
			k -= 1
			k <= 0 && return nothing

			nc = s.S[m]
			m -= 1
		end

		s.A[k] = s.S[m]
		s.S[m] = nc - 1

		if k < p - 3
			k += 1
		else
			m -= 1
		end
	end

	return nothing
end


function smith_prep!(B, C, val, a, b, c)
	if (b > p)
		val[i] += 1
		B[a,i] = c
	else
		C[:,i] .+= x[:,b] .* c
	end
end

function smith_optimize!(s, τ::Float64)
	# τ: a small positive number */
	# finds better coordinates XX[p+1..p+k1][] for the k1 Steiner points
	# of tree T by: doing a relaxation iteration. Assumes that edge lengths of old tree
	# have been per-stored in array s.L[][] */

	eqnstack = Vector{Int}(undef, n)
	leafQ    = Vector{Int}(undef, n)
	val      = Vector{Int}(undef, n)

	lqp = eqp = 0
	k1 = s.p - 2

	# First: compute B array, C array, and valences. Set up leafQ.
	for i = k1+1:-1:1
		n0 = adj[1,i]
		n1 = adj[2,i]
		n2 = adj[3,i]

		q0 = 1.0 / (s.L[1,i] + τ)
		q1 = 1.0 / (s.L[2,i] + τ)
		q2 = 1.0 / (s.L[3,i] + τ)

		# printf("q: %20.20g %20.20g %20.20g\n", q0, q1,  q2);

		# Have now figured out reciprocal distances q0, q1, q2 from
		# Steiner pt. i to neighbors n0, n1, n2 */
		t = q0 + q1 + q2

		q0 /= t
		q1 /= t
		q2 /= t

		val[i] = 0

		B[:,i] .= 0.0
		C[:,i] .= 0.0

		smith_prep!(B, C, val, 0, n0, q0)
		smith_prep!(B, C, val, 1, n1, q1)
		smith_prep!(B, C, val, 2, n2, q2)

		#printf("SP: %20.20g %20.20g\n", XX[i + p][0], XX[i + p][1]); */
		#printf("C: %20.20g %20.20g\n", C[i][0], C[i][1]); */

		# Now: Steiner point i has Steiner valence val[i];
		# coords obey eqns XX[i+p][] = sum(j)of B[j,i]*XX[nj][] + C[i][] */
		if (val[i] <= 1)
			leafQ[lqp] = i
			lqp += 1 # puts leafs on leafQ
		end
	end

	# Have set up equations - now to solve them.
	# Second: eliminate leaves
	while lqp > 1
		lqp -= 1
		i = leafQ[lqp]
		val[i] -= 1
		i2 = i + p

		# now to eliminate leaf i

		eqnstack[eqp] = i
		eqp+=1
		# push i onto stack
		for j = 1:3
			if !iszero(B[j,i])
				break # neighbor is #j
			end
		end

		q0 = B[j,i]
		j = adj[j,i] - p # neighbor is j
		val[j] -= 1

		if (val[j] == 1)
			leafQ[lqp] = j
			lqp += 1
		end # new leaf?
		
		for m = 1:3
			if (adj[m,j] == i2)
				break
			end
		end

		q1 = B[m,j]
		B[m,j] = 0.0

		t = 1.0 - q1 * q0
		t = 1.0 / t;
		
		
		B[:,j] .*= t
		C[:,j] .= t .* (C[:,j] .+ q1 .* C[:,i])
	end

	# III: Solve 1-vertex tree!
	i = leafQ[1]
	i2 = i + p

	x[:,i2] .= C[:,i]

	# /* printf("SP(1): %20.20g %20.20g\n", XX[i + p][0], XX[i + p][1]); */

	# IV: backsolve
	while eqp > 0
		eqp -= 1
		i = eqnstack[eqp]
		i2 = i + p

		for j = 1:3
			if !iszero(B[j,i])
				break # neighbor is #j
			end
		end

		q0 = B[j,i]
		j = adj[j,i] # neighbor is j */

		x[:,i2] .= C[:,i] + q0 * x[:,j]
	end

	# /* printf("SP(2): %20.20g %20.20g\n", XX[i + p][0], XX[i + p][1]); */

	return
end

function smith_build_tree!(s::SmithData{T,U}, k::Integer) where {T,U}
    # builds tree represented by topvec[1..k]. Initial location of new Steiner pts
    # is a slightly random perturbation of the centroid of its neighbors
	
    # First build the tree corresponding to the null vector
	m = s.p + 1

	s.ℵ[1:3,1] .= 1:3
	s.E[1,1:3] .= 1:3
	s.E[2,1:3] .= m

	for i = 1:s.N
		s.x[i,m] = sum(s.x[i,1:3]) / 3.0 + ε₂ * s.α * rand(s.rng)
	end

	for i = 2:k+1
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
			s.x[j,sn] = (s.x[j,ea] + s.x[j,eb] + x[j,en]) / 3.0 + ε₂ * α * rand(rng)
		end
	end

	s.k = k + 3
	
	# Tree is now built. Initial coords in general position.
	
	return
end

function smith_dist(s::SmithData{T,U}, i::Integer, j::Integer) where {T,U}
	return norm(s.x[:, i] - s.x[:, j])
end

function smith_length_delta!(s::SmithData{T,U}, m::Integer, i::Integer, k::Integer) where {T,U}
	if u[m] < k
		Δℓ = smith_dist(s, s.u[m], k)
		
		s.L[m,i] = Δℓ
		s.u[m] -= s.p

		if s.u[m] > 0
			for j = 1:3
				if s.ℵ[j,s.u[m]] == k
					s.L[j,s.u[m]] = Δℓ
					break
				end
			end
		end

		return Δℓ
	else
		return zero(T)
	end
end

function smith_length!(s::SmithData{T,U}) where {T,U}
	# Stores edge lengths of tree T in array s.L[1..k1][0..2] and returns total length.
	m = s.k - 2
	ℓ = zero(T)

	for i = 1:m
		s.u[:] .= s.ℵ[:, i]

		ℓ += smith_length_delta!(s, m, i, i + p)
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
		s.v[:,:] .= s.x[m,s.u] .- s.x[m,k]

		for (l₁, l₂) in [(1,2), (1,3), (2, 3)]
			let Δε = 2 * (v[:,l₁]' * v[:,l₂])
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
