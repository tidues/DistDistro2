import randist as rt

# read data
gname = 'manhattan'
fpath = '../data/'
g = rt.readGraph(fpath, gname)
rt.get_largest_component(g)
