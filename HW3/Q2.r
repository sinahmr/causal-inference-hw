library(igraph)

remove_directed = function(adj) {
        p = length(adj[1,])
        for (i in 1:p) {
                for (j in 1:p) {
                        if (adj[i, j] == 0) {
                                adj[j, i] = 0
                        }
                }
        }
        return(adj)
}

get_comps = function(adj) {
        g  = graph.adjacency(adj)
        clu = components(g)
        return(groups(clu))
}

directed_chain = function(adj, v) {
        if (length(adj) == 1) return(1)
        
        p = length(adj[1,])
        a = list(v)
        b = (1:p)[-v]
        
        while(length(b) != 0) {
                t = array(0, dim=p)
                for(i in a) {
                        for(j in b) {
                                if (adj[i, j] == 1 && adj[j, i] == 1) {
                                        t[j] = 1
                                        adj[j, i] = 0
                                }
                        }
                }
                
                t = which(t == 1)
                flag = 1
                while(flag == 1) {
                        flag = 0
                        for (x in 1:p) {
                                for (y in t) {
                                        for (z in t) {
                                                if (adj[x, y] == 1 && adj[y, x] == 0 && adj[y, z] == 1 && adj[z, y] == 1 && adj[x, z] == 0 && adj[z, x] == 0) {
                                                        adj[z, y] = 0
                                                        flag = 1
                                                }
                                        }
                                }
                        }
                }
                
                a = t
                b = setdiff(b, t)
        }
        return(adj)
}

size_mec = function(adj) {
        if (length(adj) == 1) return(1)
        
        p = length(adj[1,])
        n = sum(adj) / 2
        
        if (n == p - 1) return(p)
        if (n == p) return(2 * p)
        if (n == p*(p-1)/2 - 2) return((p*p - p - 4) * factorial(p - 3))
        if (n == p*(p-1)/2 - 1) return(2 * factorial(p - 1) - factorial(p - 2))
        if (n == p*(p-1)/2) return(factorial(p))
        
        s = 0
        for (j in 1:p) {
                chain_adj = directed_chain(adj, j)
                chain_adj = remove_directed(chain_adj)
                comps = get_comps(chain_adj)
                pi = 1
                for (comp in comps) {
                        pi = pi * size_mec(chain_adj[comp, comp])
                }
                s = s + pi
        }
        return(s)
}

main = function(adj) {
        adj = remove_directed(adj)
        comps = get_comps(adj)
        pi = 1
        for (comp in comps) {
                sub_adj = adj[comp, comp]
                pi = pi * size_mec(sub_adj)
        }
        return(pi)
}


# result: 8
toy_adj_1 = rbind(
        c(0, 0, 0, 0, 0),
        c(1, 0, 1, 1, 0),
        c(1, 1, 0, 1, 0),
        c(0, 1, 1, 0, 1),
        c(0, 0, 0, 1, 0)
)
print(main(toy_adj_1))

# result: 4
toy_adj_2 = rbind(
        c(0, 1, 0, 1, 0),
        c(1, 0, 1, 0, 0),
        c(0, 0, 0, 0, 1),
        c(1, 0, 0, 0, 1),
        c(0, 0, 0, 1, 0)
)
print(main(toy_adj_2))

# result: 30
toy_adj_3 = rbind(
        c(0, 1, 0, 0, 0),
        c(1, 0, 1, 1, 1),
        c(1, 1, 0, 1, 1),
        c(1, 1, 1, 0, 1),
        c(1, 1, 1, 1, 0)
)
print(main(toy_adj_3))
