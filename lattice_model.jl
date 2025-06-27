#=
3D lattice gas model, in the presence of an Ising surface (membrane), coupled via tethers.


Hamiltonian:

H = -Jb*sum(ni*nj) -Jm*sum(si*sj)-Jt*sum(ni)-mu*sum(ni),

where the first and fourth terms are for bulk.
The second term sums over all the Ising spins in the membrane.
The third term is summed over all the indices where tethers exist.

Set of possible moves:
1. Particle displacement (swap)
2. Particle addition/deletion
3. Membrane lipid movement (swap)
4. Tether movement
=#

using Distributions, JLD, StatsBase, Random, Base.Threads, Plots


# Model parameters
const L = 20		# box size
const Jb = 0.5     	# bulk coupling (kT, Jbc~0.88625)
const Jm = 2.0          # membrane spins coupling constant (Kt, Jmc~0.441)
const Jt = 5.0          # coupling b/w tether and particle 
const mu = -1.6		# bulk chemical potential
const Td = 0.5          # tether density (fraction of up-spins having tethers)
const Lt = 5            # tether length
const phi_mem =  0.5    # membrane lipid composition

# Implementation parameters
const Nsweeps = 5*10^4
const Nsteps = (L^3 + L^2)      # number of MC steps per sweep


"""
Returns an array of particles; initialized with random values from {0,1}
"""
function bulk_initial(l::Int64, phi=0.5::Float64)::Array{Int8}
    weights = [1-phi, phi]
    temp = sample([Int8(0), Int8(1)], Weights(weights), (l,l,l))
    return temp
end

"""
Initializes (randomly) Ising spins on the membrane. Values from {-1, 1}
"""
function mem_initial(l::Int64, phi=0.5::Float64)::Matrix{Int8}
    temp = -1*ones(l,l)
    # randomly place +1 sites on phi fraction of sites
    nup = Int64(round(phi*l*l))
    up_sites = shuffle!([[n,m] for n in 1:l, m in 1:l])[1:nup]
    for n in up_sites
        temp[n[1], n[2]] = Int8(1)    
    end
    return temp
end

"""
Initializes tethers. Tether density Td is the fraction of up spins having a tether attached. 
"""
function teth_initialize(memb::Matrix{Int8})::Vector{Vector{Int8}}
    l = size(memb)[1]
    up_pos = findall(e-> e ==1, memb)                   # positions of all the +1 spins on the membrane
    num_teth = Int8(round(Td*length(up_pos)))                 # number of tethers
    selected_upspin_pos = shuffle!(up_pos)[1:num_teth]        # randomly selects positions to place tethers
    temp = [Int8.([n[1], n[2]]) for n in selected_upspin_pos]
    return temp
end

"""
Returns list of neighbours of a site (for bulk)
"""
function neighbours(ind::Vector{Int64}, l::Int64)
	neigh = Vector{Vector{Int64}}(undef, 6)
	i,j,k = ind
	# PBC for x direction exists
	neigh[1] = [mod1(i+1, l), j, k]
	neigh[2] = [mod1(i-1, l), j, k]
	# PBC for y direction exists
	neigh[3] = [i, mod1(j+1, l), k]
	neigh[4] = [i, mod1(j-1, l), k]
	# no PBC for z direction
    if k == l
		neigh[5] = [i,j,l-1]
		return neigh[1:5]
	elseif k == 1
		neigh[5] = [i,j,2]
		return neigh[1:5]
	elseif ((k>=2) && (k<=l-1))
		neigh[5] = [i, j, k+1]
		neigh[6] = [i, j, k-1]
		return neigh
	end
end
	

""" 
Calculates the total energy of the system. All three components.
"""
function calculate_energy(bulk::Array{Int8}, membrane::Matrix{Int8}, tethers::Vector{Vector{Int8}})::Float64
	#bulk
    l = size(bulk)[1]
    H1_matrix = zeros(Float64, size(bulk))
	for k in 1:l
		for j in 1:l
			for i in 1:l
                if bulk[i,j,k] == 0
                    H1_matrix[i,j,k] = 0
                else
                    H1_matrix[i,j,k] = -0.5*Jb*sum([bulk[n[1], n[2], n[3]] for n in neighbours([i,j,k], l)])
                    H1_matrix[i,j,k] += -Jt*([i,j] in tethers)*(k>=l-Lt)
                end
			end
		end
	end
	bulk_energy = round((sum(H1_matrix) -mu*sum(bulk)), digits=6)

    #membrane
    H2_matrix = zeros(Float64, size(membrane))
    for j in 1:l
        for i in 1:l
            neigh_pts = [[mod1(i+1,l), j], [mod1(i-1,l), j], [i, mod1(j+1,l)], [i, mod1(j-1,l)]]
            H2_matrix[i,j] = -0.5*Jm*membrane[i,j]*sum([membrane[n[1], n[2]] for n in neigh_pts])
        end
    end
    memb_energy = sum(H2_matrix)
    total_energy = round(bulk_energy+memb_energy, digits=6)
    return total_energy#, H1_matrix, H2_matrix
end




##################### Possible moves ####################################################
"""
Does a particle displacement (global, since we are interested in eql. behaviour)
"""
function displace!(bulk::Array{Int8}, tethers::Vector{Vector{Int8}})::Tuple{Array{Int8, 3}, Float64}
	l = size(bulk)[1]
	k = rand(1:l, 3)          # choose random sites for global exchange            
	q = rand(1:l, 3)
	dE = 0.0
	# making the step
	if bulk[k[1], k[2], k[3]] != bulk[q[1], q[2], q[3]]
		kn = filter(e-> e != q, neighbours(k,l))				#neighbours of k, except q
		qn = filter(e-> e != k, neighbours(q,l))			#neighbours of q, except k
        # for site k
		Hk = -Jb*sum([bulk[n[1], n[2], n[3]] for n in kn])
		Hk += -Jt*([k[1], k[2]] in tethers)*(k[3]>= l-Lt)       #contribution from tether
        
        # for site q
		Hq = -Jb*sum([bulk[n[1], n[2], n[3]] for n in qn])
		Hq += -Jt*([q[1], q[2]] in tethers)*(q[3]>= l-Lt)       #contribution from tether
		
        #total change in energy due to displacement
        dE = round((Hk-Hq)*(-1)^(bulk[k[1], k[2], k[3]]), digits=6)     
				
        if rand() < minimum([1, exp(-dE)])	# acceptance prob.
			bulk[k[1], k[2], k[3]] = Int8(mod(bulk[k[1], k[2], k[3]]+1, 2))
			bulk[q[1], q[2], q[3]] = Int8(mod(bulk[q[1], q[2], q[3]]+1, 2))
			return bulk, dE
		else								# if Monte Carlo step is rejected
			return bulk, 0.0
	    end
	else									# if both sites have the same value
		return bulk, 0.0
	end
end


"""
Adds or removes a particle
"""
function particle!(bulk::Array{Int8}, tethers::Vector{Vector{Int8}})::Tuple{Array{Int8, 3}, Float64}
	l = size(bulk)[1]
	k = rand(1:l, 3)          # choose random site in bulk              
	if rand() < 0.5								
    # particle addition
        if bulk[k[1], k[2], k[3]] == 1
            return bulk, 0.0
        else
            dE = -Jb*sum([bulk[n[1], n[2], n[3]] for n in neighbours(k, l)])    # nearest neighbour contribution
            dE += -Jt*([k[1], k[2]] in tethers)*(k[3]>= l-Lt) - mu              # contribution from tether and chempot
            accprob = round(exp(-dE), digits=6)
            if rand() < minimum([1,accprob])
                bulk[k[1], k[2], k[3]] = Int8(1)
                return bulk, dE
            else
                return bulk, 0.0
            end
        end
	else 
    # particle deletion
        if bulk[k[1], k[2], k[3]] == 0
            return bulk, 0.0
        else
            dE = Jb*sum([bulk[n[1], n[2], n[3]] for n in neighbours(k, l)])    # nearest neighbour contribution
            dE += Jt*([k[1], k[2]] in tethers)*(k[3]>= l-Lt) + mu              # contribution from tether and chempot
            accprob = round(exp(-dE), digits=6)
            if rand() < minimum([1,accprob])
                bulk[k[1], k[2], k[3]] = Int8(0)
                return bulk, dE
            else
                return bulk, 0.0
            end
        end
    end
end


"""
Global exchange b/w two membrane lipids (only those without tethers)
"""
function lipid_move!(memb::Matrix{Int8}, tethers::Vector{Vector{Int8}})
    l = size(memb)[1]
    # select two site not occupied by tethers
    k = rand(1:l, 2)
    while k in tethers
        k = rand(1:l, 2)
    end
    q = rand(1:l, 2)
    while q in tethers
        q = rand(1:l, 2)
    end

    if memb[k[1], k[2]] != memb[q[1], q[2]]
        kn = filter(e-> e != q, [[mod1(k[1]+1,l),k[2]], [mod1(k[1]-1,l),k[2]], [k[1],mod1(k[2]+1,l)], [k[1], mod1(k[2]-1,l)]])			# neighbours of k, except q
		qn = filter(e-> e != k, [[mod1(q[1]+1,l),q[2]], [mod1(q[1]-1,l),q[2]], [q[1],mod1(q[2]+1,l)], [q[1], mod1(q[2]-1,l)]])			# neighbours of q, except k
        dE = 2*Jm*( memb[k[1], k[2]]* sum([memb[n[1], n[2]] for n in  kn]) + memb[q[1], q[2]]*sum([memb[n[1], n[2]] for n in qn]))
        accprob = round(exp(-dE), digits=6)
        # Metropolis step
        if rand() < minimum([1,accprob])
            memb[k[1], k[2]] *= -1
            memb[q[1], q[2]] *= -1
            return memb, dE
        else
            return memb, 0.0
        end
    else
        return memb, 0.0
    end
end


"""
Moves tether to any neighbouring site with +1 states. 
"""
function tether_move!(bulk::Array{Int8}, memb::Matrix{Int8}, tethers::Vector{Vector{Int8}})
    l = size(memb)[1]
    k = rand(tethers)               # select tether to move
    neighbouring_lipids = [[mod1(k[1]+1,l),k[2]], [mod1(k[1]-1,l),k[2]], [k[1],mod1(k[2]+1,l)], [k[1], mod1(k[2]-1,l)]]
    possible_sites = filter(e-> memb[e[1], e[2]] == 1, neighbouring_lipids)        # filtering neighbours with +1 state
    possible_sites = filter!(e-> !(e in tethers), possible_sites)
    if length(possible_sites) == 0
        return tethers, 0.0
    else
        q = rand(possible_sites)                # neighbouring site to move tether on to
        dE = Jt*(sum([bulk[k[1], k[2], z] for z in l-Lt:l]) - sum([bulk[q[1], q[2], z] for z in l-Lt:l]))
        accprob = round(exp(-dE), digits=6)
        # Metropolis step
        if rand() < minimum([1,accprob])
            tethers = filter!(e-> e != k, tethers)
            tethers = push!(tethers, q)
            return tethers, dE
        else
            return tethers, 0.0
        end
    end
end

############################ Execution ##################################################

""" 
Carries out a Markov step 
"""
function markov_step!(bulk::Array{Int8}, memb::Matrix{Int8}, tethers::Vector{Vector{Int8}})
    weights = [30, 30, 20, 20]                # Weights for choosing particle displacement, add/delete, tether move and lipid move
    # choosing random step to carry out
    process = sample([1,2,3,4], Weights(weights))
    if process == 1             # particle displacement
        bulk, dE = displace!(bulk, tethers) 
    elseif process == 2         # particle addition/deletion
        bulk, dE = particle!(bulk, tethers)
    elseif process == 3         # tether movement
        tethers, dE = tether_move!(bulk, memb, tethers)
    elseif process == 4         # membrane lipid movement
        memb, dE = lipid_move!(memb, tethers)
    end
	return bulk, memb, tethers, dE
end


""" 
Main execution loop 
"""
function main()
    #initializing
    bulk = bulk_initial(L, 0.5)
    memb = mem_initial(L, phi_mem)
    tethers = teth_initialize(memb)

    # storing the states (bulk + membrane + tethers + energy value) after each sweep
	energy = zeros(Nsweeps+1)
	bulk_data = Array{Int8}(undef, (L,L,L,Nsweeps+1))
    membrane_data = Array{Int8}(undef, (L, L, Nsweeps+1))
    tether_data = Vector{Vector{Vector{Int8}}}(undef, Nsweeps+1)

    # storing the initial values
	en = calculate_energy(bulk, memb, tethers)[1]
	energy[1] = en
	bulk_data[:,:,:,1] = bulk
    membrane_data[:,:,1] = memb
    tether_data[1] = tethers

    # Monte Carlo steps
	for n in 2:Nsweeps+1
		for t in 1:Nsteps
            bulk, memb, tethers, dE = markov_step!(bulk, memb, tethers)
            en += dE
        end
        energy[n] = en
        bulk_data[:,:,:,n] = bulk
        membrane_data[:,:,n] = memb
        tether_data[n] = tethers
	end
	return bulk_data, energy, membrane_data, tether_data
end
