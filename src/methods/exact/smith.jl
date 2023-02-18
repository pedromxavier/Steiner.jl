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

function _length()
	# Stores edge lengths of tree T in array EL[1..k1][0..2] and returns total length. */
	
	leng = 0.0
	k1 = N - 2
	
	for i = 2:k1+1
		i2 = i + p

		n0 = adj[1,i]
		n1 = adj[2,i]
		n2 = adj[3,i]

		if n0 < i2
			t = norm(x[:,n0], x[:,i2])
			leng += t
			EL[1,i] = t
			n0 -= p

			if n0 > 0
				for j = 1:3
					if adj[j,n0] == i2
						EL[j,n0] = t
						break
					end
				end
			end
		end

		if n1 < i2
			t = norm(x[:,n1], x[:,i2])
			leng += t
			EL[2,i] = t
			n1 -= p
			
			if n1 > 0
				for j = 1:3
					if adj[j,n1] == i2
						EL[j,n1] = t
						break
					end
				end
			end
		end

		if n2 < i2
			t = norm(x[:,n2], x[:,i2])
			leng += t
			EL[3,i] = t
			n2 -= p
			
			if n2 > 0
				for j = 1:3
					if adj[j,n2] == i2
						EL[j,n2] = t
						break
					end
				end
			end
		end
	end 

	# Have now figured out distance EL[i][00.3] from Steiner pt. i to neighbors. */

	return leng
end

function _error()
	# Returns the error figure of tree T with Steiner coords in XX[][].
	# Assumes edge lengths have been pre-stored in array EL[][]. */
	
	k1 = N - 2
	efig = 0.0

	for i = 2:k1+1
		i2 = i + p
		n0 = adj[1,i]
		n1 = adj[2,i]
		n2 = adj[3,i]

		d12 = d01 = d02 = 0.0

		for m = 1:N
			t = x[m,i2]
			r = x[m,n0] - t
			s = x[m,n1] - t
			t = x[m,n2] - t

			d12 += s*t; d01 += r*s; d02 += r*t;
		end

		# only angles < 120 cause error */

		t = d12 + d12 + EL[2,i] * EL[3,i]

		if (t > 0.0)
			efig += t
		end

		t = d01 + d01 + EL[1,i] * EL[2,i]

		if (t > 0.0)
			efig += t
		end

		t = d02 + d02 + EL[1,i] * EL[3,i]

		if (t > 0.0)
			efig += t
		end

	end
	
	return sqrt(efig)
end

function solve(::SmithMethod, points::Matrix{K}) where {K}
    N = size(points, 1) # space dimension
    p = size(points, 2) # terminal points
	s = p - 2
	n = p + s

	# Inputs p, N, sites; outputs succesive best Steiner
    # trees found. Best tree's topology-vector is stored in BESTVEC.	
	A = Vector{Int}(undef, n)
	x = Matrix{Float64}(undef, N, 2n)
	
	rng = Random.Xoshiro()

	# (1) Input

	for i = 2:p+1, j = 1:N
		x[j,i] = ...
	end

	# Constants
	α = 0.0     # Scale
	β = 0.99999
	δ = 0.0001
	

	for i = 2:p+1, j = 1:N
		q = x[j,i] - x[j,2]

		if q < 0.0
			q = -q
		end

		if q > α
			α = q
		end
	end

	if p == 3 # Deal with special case of 3 sites
		_buildtree(0, A)

		q = _length()
		r = _error()

		while true
			optimize(δ * r / p)
			q = _length()
			r = _error()

			if !(r > q * δ)
				break
			end
		end

		output_tree()

		return
	end
	
	# (2) Preprocessing and initialization
	# Optionally, sort sites in some nice order here
	# ℓ = any upper bound ont the length of the SMT
	ℓ = HUGE
	k    = 1
	m    = 0
	ct   = 0

	# ct counts backtrack iters. Unused at present
	while true
		# 3: candidate leaf generation and backtracking */
		nc = 0
		ct += 1

		# if (ct % 10000 == 0) { printf("still running (%ld)!\n",ct); fflush(stdout); }
		
		l = 2 * k + 1

		while l > 0 && k <= p - 3
			# Build the tree represented by the topol. vector A[1..k]
			A[k] = l

			_buildtree(k, A)

			# .. and optimize it until either obviously bad or small error figure happens
			@label ITER

			q = _length()
			r = _error()

			if (q - r < ℓ)
				if (r > 0.005 * q)
					optimize(δ * r / p)
					@goto ITER
				end

				if (k >= p-3)
					while true
						optimize(δ * r / p)
						q = _length()
						r = _error()

						if !(r > q* δ)
							break
						end
					end

					if (q < ℓ)
						for i = 2:k+1
							BESTVEC[i] = A[i]
						end

						if (q < ℓ * β)
							output_tree()
						end

						ℓ = q;
					end
				else
					i = nc
					nc += 1

					while (i > 1 && LEN[i] < q)
						STACK[m+i+1] = STACK[m+i]
						LEN[i+1] = LEN[i]
						i -= 1
					end

					i += 1

					STACK[m+i] = l
					
					LEN[i] = q
				end
			end

			l -= 1
		end

		m = m + nc;

		while (nc <= 0)
			k -= 1
			if (k <= 0)
				return 0 # exit(0)
			end

			nc = STACK[m]
			m -= 1
		end

		A[k] = STACK[m]
		STACK[m] = nc - 1
		
		if (k < p - 3)
			k += 1
		else
			m -= 1
		end
	end

	return 0
end

function _buildtree(k::Int, topvec::Vector{Int})
	rng = Random.Xoshiro()

	# Tiny factor
	ϵ = 0.001

	# Scale
	α = 1.0

    # builds tree represented by topvec[1..k]. Initial location of new Steiner pts
    # is a slightly random perturbation of the centroid of its neighbors
	
    # First build the tree corresponding to the null vector
    N = 3
	m = p+1

	for i = 1:3
		adj[i,1]  = i
		edge[1,i] = i
		edge[2,i] = m
	end

	for i = 1:N
		XX[i,m] = (XX[i,1] + XX[i,2] + XX[i,3]) / 3.0 + ϵ * α * rand(rng)
	end

	for i = 2:k+1
		# Now: do vector element topvec[i]
		en = i + 3
		m  = i + 1
		sn = m + p
		e  = topvec[i]
		ea = edge[1,e]
		eb = edge[2,e]
		adj[1,m] = ea
		adj[2,m] = eb
		adj[3,m] = en
		m = ea - p

		if m > 0
			for j = 1:3
				if adj[j,m] == eb
					adj[j,m] = sn
					break
				end
			end
		end

		m = eb - p

		if m > 0
			for j = 1:3
				if adj[j,m] == ea
					adj[j,m] = sn
					break
				end
			end
		end

		edge[2,e] = sn
		
		e = 2en - 4

		edge[1,e] = en
		edge[2,e] = sn

		e += 1

		edge[1,e] = eb
		edge[2,e] = sn

		for j = 1:N
			x[j,sn] = (x[j,ea] + x[j,eb] + x[j,en]) / 3.0 + ϵ * α * rand(rng)
		end
	end

	N = k + 3
	
	# Tree is now built. Initial coords in general position.
	
	return
end

"""
/* Global variables: */
double ℓ, α, N;
int p, N;
static int BESTVEC[N], STACK[N*N], adj[N-2][3], edge[2*N][2];
static double XX[N*2][MAXDIM], LEN[N], EL[N][3];
"""

function _prep!(B, C, val, a, b, c)
	if (b > p)
		val[i] += 1
		B[a,i] = c
	else
		C[:,i] .+= x[:,b] .* c
	end
end

function _optimize(τ::Float64)
	# τ: a small positive number */
	# finds better coordinates XX[p+1..p+k1][] for the k1 Steiner points
	# of tree T by: doing a relaxation iteration. Assumes that edge lengths of old tree
	# have been per-stored in array EL[][] */

	B = Matrix{Float64}(undef, 3, n)
	C = Matrix{Float64}(undef, N, n)

	eqnstack = Vector{Int}(undef, n)
	leafQ    = Vector{Int}(undef, n)
	val      = Vector{Int}(undef, n)

	lqp = eqp = 0
	k1 = p - 2

	# First: compute B array, C array, and valences. Set up leafQ.
	for i = k1+1:-1:1
		n0 = adj[1,i]
		n1 = adj[2,i]
		n2 = adj[3,i]

		q0 = 1.0 / (EL[1,i] + τ)
		q1 = 1.0 / (EL[2,i] + τ)
		q2 = 1.0 / (EL[3,i] + τ)

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

		_prep!(B, C, val, 0, n0, q0)
		_prep!(B, C, val, 1, n1, q1)
		_prep!(B, C, val, 2, n2, q2)

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
	while(lqp > 1)
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