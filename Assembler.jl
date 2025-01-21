using BioSequences
using LightGraphs
using GraphPlot
using Base.Threads

# construct the De-Bruijn Graph
function construct_de_bruijn_graph(sequences::Vector{String}, k::Int)
    graph = Dict{String, Vector{String}}()
    for sequence in sequences
        # generate k-mers for each sequence 
        for i in 1:(length(sequence) - k + 1)
            prefix = sequence[i:i+k-2]
            suffix = sequence[i+1:i+k-1]
            if !haskey(graph, prefix)
                graph[prefix] = []
            end
            push!(graph[prefix], suffix)
        end
    end
    return graph
end

# visualize graph in adjacency list format
function print_de_bruijn_graph(graph::Dict{String, Vector{String}})
    for (node, edges) in graph
        println("$node -> ", join(edges, ", "))
    end
end

function find_start_kmer(debruijn_graph::Dict{String, Vector{String}})
    # look for kmer with no predecessor
    for (kmer, neighbors) in debruijn_graph
        if !any(kmer in successors for successors in values(debruijn_graph))
            return kmer 
        end
    end

    # no kmer with no predecessor, find one with fewer predecessors than successors
    # hold degrees
    in_degree = Dict{String, Int}()
    out_degree = Dict{String, Int}()

    # init the degrees
    for kmer in keys(debruijn_graph)
        in_degree[kmer] = 0
        out_degree[kmer] = length(debruijn_graph[kmer])
    end

    # iterate through edges and find degrees
    for neighbors in values(debruijn_graph)
        for neighbor in neighbors
            in_degree[neighbor] = get(in_degree, neighbor, 0) + 1
        end
    end

    # find k-mer with fewer predecessors than successors
    for kmer in keys(debruijn_graph)
        if in_degree[kmer] < out_degree[kmer]
            println("Start k-mer: $kmer (in-degree: $in_degree[kmer], out-degree: $out_degree[kmer])")
            return kmer 
        end
    end

    # no kmer exists, find max in-degree
    start_kmer = argmax(in_degree)
    max_in_degree = in_degree[start_kmer]
    println("Fallback Start k-mer: $start_kmer (in-degree: $max_in_degree)")
    return start_kmer
end

function find_eulerian_path(graph::Dict{String, Vector{String}}, start::String)
    path = []
    stack = [start] # use the start node from arguments
    local_graph = deepcopy(graph)  # copy to not modify

    # path search
    while !isempty(stack)
        current_node = stack[end]
        if haskey(local_graph, current_node) && !isempty(local_graph[current_node])
            # move to next and pop
            next_node = pop!(local_graph[current_node])
            push!(stack, next_node)
        else
            # backtrack and add current to the path
            push!(path, pop!(stack))
        end
    end

    return reverse(path)
end

# reconstruct sequence
function reconstruct_sequence(path::Vector{Any}, k::Int)
    # go through the path and append a letter for each node
    sequence = path[1]
    for node in path[2:end]
        sequence *= last(node)
    end
    return sequence
end

# find the shortest string length of vectors
function shortest_string_length(sequences::Vector{String})
    return minimum(map(length, sequences))
end

# run sequence for each k value
function parallelProcess(sequence::Vector{String}, k_min::Int, k_max::Int)
    results = Vector{String}(undef, k_max - k_min + 1)

    @threads for k in k_min:k_max
        println("Processing for k = $k...")    
        try
            graph = construct_de_bruijn_graph(sequence, k)
            print_de_bruijn_graph(graph)
            start_kmer = find_start_kmer(graph)
            path = find_eulerian_path(graph, start_kmer)
            reconstructed_sequence = reconstruct_sequence(path, k)
            results[k - k_min + 1] = reconstructed_sequence
        catch e
            println("Error for k = $k: $e")
            results[k - k_min + 1] = "Error"
        end
    end

    return results
end

# sequence = "AGCTCTAGCAGAGC"


# Test 1
sequences = ["GGCGTCTATATCT", "GGCGTCGATATCT", "TCTATATCTCGGCTCTAGG"]

k_min = 10
k_max = 10

# run computation for each k value
results = parallelProcess(sequences, k_min, k_max)
# graph = construct_de_bruijn_graph(sequences, k_min)
# print_de_bruijn_graph(graph)
# display results for all sequences
for k in k_min:k_max
    println("k = $k, Reconstructed Sequence: $(results[k - k_min + 1])")
end



# Test 2
sequences = ["GGCGTCTATATCT"]

k_min = 3
k_max = 3

# run computation for each k value
results = parallelProcess(sequences, k_min, k_max)

# display results for all sequences
for k in k_min:k_max
    println("k = $k, Reconstructed Sequence: $(results[k - k_min + 1])")
end



# Test 3: Inaccurate
sequences = ["GGCGTCTATATCT", "GGCGTCGATATCT", "TCTATATCTCGGCTCTAGG", "TATCTCGACTCTAGGCC", "TATCTCGACTCTAGGCCCTCA", "CTCGGCTCTAGCCCCTCATTTT", "GGCTCTAGGCCCTCATTTTTT", "CTCTAGGCCCTCAATTTTT"]

k_min = 3
k_max = shortest_string_length(sequences)


# run computation for each k value
results = parallelProcess(sequences, k_min, k_max)

# display results for all sequences
for k in k_min:k_max
    println("k = $k, Reconstructed Sequence: $(results[k - k_min + 1])")
end

# graph = construct_de_bruijn_graph(sequence, k)

# println("De Bruijn Graph:")
# print_de_bruijn_graph(graph)


# start_kmer = find_start_kmer(graph)
# println("Start k-mer: $start_kmer")

# path = find_eulerian_path(graph, start_kmer)
# print(path)
# reconstructed_sequence = reconstruct_sequence(path, k)

# println("Original Sequence: $sequence")
# println("Reconstructed Sequence: $reconstructed_sequence")




#-------------------------Unused Functions------------------------------
# function construct_de_bruijn_graph_multi(sequences::Vector{String}, k::Int)
#     graph = Dict{String, Vector{String}}()
#     for sequence in sequences
#         # Ensure the sequence is long enough for the given k
#         if length(sequence) < k
#             println("Sequence too short for k-mer size: $sequence")
#             continue
#         end
#         for i in 1:(length(sequence) - k)
#             kmer1 = sequence[i:i+k-1]
#             kmer2 = sequence[i+1:i+k]
#             if haskey(graph, kmer1)
#                 push!(graph[kmer1], kmer2)
#             else
#                 graph[kmer1] = [kmer2]
#             end
#         end
#     end
#     return graph
# end

# Merged into Create Graph
# generate k-mers from each sequence 
# function generate_kmers(sequence::String, k::Int)
#     kmers = []
#     for i in 1:(length(sequence) - k + 1)   # 
#         push!(kmers, sequence[i:i+k-1])
#     end
#     return kmers
# end

# Originally used to ensure a complete graph with connections to everything
# function ensure_complete_graph(graph::Dict{String, Vector{String}})
#     # Collect all unique nodes from keys and values
#     all_nodes = Set{String}(keys(graph))
#     for values in values(graph)
#         append!(all_nodes, values)
#     end

#     # Ensure every node exists as a key
#     for node in all_nodes
#         if !haskey(graph, node)
#             graph[node] = String[]  # Initialize with no outgoing edges
#         end
#     end
# end

# function reconstruct_sequence_multi(paths::Vector{Vector{Int}}, k::Int)
#     sequence = ""
#     for path in paths
#         # Decode and combine paths
#         sub_sequence = reconstruct_sequence(path, k)
#         sequence *= sub_sequence[k:end]  # Append avoiding duplicate overlap
#     end
#     return sequence
# end
