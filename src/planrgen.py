
def gen_planar(side_nodes = 10):
    text = ''
    maxv = side_nodes ** 2
    vers = range(1, maxv + 1)
    text = 'i\tj\tl\tx\ty\n'
    for v in vers:
        for u in nbr(v, side_nodes):
            text += (
                    str(v) + '\t' + 
                    str(u) + '\t' + 
                    str(1) + '\t' + 
                    str(1.0/maxv) + '\t'+ 
                    str(1.0/maxv) + '\n'
                    )

    wf = open('planar_side_'+str(side_nodes)+'.dat', 'w')
    wf.write(text)

# find neigbors given a node
def nbr(v, side_nodes):
    maxv = side_nodes ** 2
    nbrs = []

    if v % side_nodes != 0:
        nbrs.append(v+1)

    if v + side_nodes <= maxv:
        nbrs.append(v+side_nodes)

    return nbrs


gen_planar(side_nodes=5)


