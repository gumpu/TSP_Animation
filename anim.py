#!/usr/bin/python2

# This only runs under Python2 

import math
import time
import random
import datetime
from node import Node
from node import length
from node import total_length
from node import dump
import matplotlib.pyplot as plt
import numpy

l_min = None
pic = 0
# v = [-100,4100,-100,2100]

def framec( solution, nc ):
    global pic
    cities = [ (n.x,n.y) for n in solution ]
    cities = numpy.array( cities )
    plt.axis( [-100,4100,-100,2100] )
    plt.axis('off')
    plt.plot(cities[:,0],cities[:,1],'ko')
    plt.title('{} Cities, 5 Search Algorithms'.format(nc))
    plt.savefig( ("%05d" % pic)+'.png')
    plt.clf()
    print pic
    pic += 1

def frame0(solution, nodes, l, title):
    global pic
    cities = [(n.x,n.y) for n in solution]
    cities.append(( solution[0].x, solution[0].y ))
    cities = numpy.array( cities )
    all_node = [(n.x,n.y) for n in nodes]
    all_nodes = numpy.array(all_node)

    plt.axis([-100, 4100, -100, 2100])
    plt.axis('off')
    plt.plot(cities[:,0], cities[:,1], 'bo-')
    plt.plot(all_nodes[:,0], all_nodes[:,1], 'ko')
    plt.title('{}  Tour length {:.1f}'.format(title,l))
    plt.savefig( ("%05d" % pic)+'.png')
    plt.clf()
    print pic
    pic += 1

nn = 0
def frame(nodes, solution, sn, t, c, y, x, z, gain):
    global pic
    global nn

    cities = [(n.x,n.y) for n in solution]
    cities = numpy.array(cities)

    cities2 = [(c.x,c.y), (y.x,y.y)]
    cities3 = [(x.x,x.y), (z.x,z.y)]
    cities2 = numpy.array(cities2)
    cities3 = numpy.array(cities3)

    plt.plot(cities[:,0],cities[:,1],'bo-')
    #plt.scatter(cities[:,0], cities[:,1],s=50,c='k')

    if gain < 0:
        plt.scatter(cities2[:,0], cities2[:,1], c='r',s=180)
        plt.plot(cities2[:,0],cities2[:,1], c='r',linewidth=2)
        plt.scatter(cities3[:,0], cities3[:,1],c='b',s=150)
        plt.plot(cities3[:,0],cities3[:,1], c='r',linewidth=2)

    else:
        plt.scatter(cities2[:,0], cities2[:,1], c='g',s=180)
        plt.plot(cities2[:,0], cities2[:,1], c='g',linewidth=2)
        plt.scatter(cities3[:,0], cities3[:,1], c='b',s=150)
        plt.plot(cities3[:,0], cities3[:,1], c='g',linewidth=2)

    plt.axis( [-100,4100,-100,2100] )
    plt.axis('off')

    plt.title('(4)  SA Temp {:4.1f} Best Tour {:6.1f}\nSwaps {}  Gain {:12.2f} '.format(t,l_min,nn,gain))
    plt.savefig( ("%05d" % pic)+'.png')
    plt.clf()
    pic += 1
    print pic


def frame4(nodes, solution, sn, c, y, x, z, gain):
    global pic
    global nn

    l_min = total_length( nodes, solution )
    cities = [ (n.x,n.y) for n in solution ]
    cities = numpy.array( cities )

    cities2 = [ (c.x,c.y), (y.x,y.y) ]
    cities3 = [ (x.x,x.y), (z.x,z.y) ]
    cities2 = numpy.array( cities2 )
    cities3 = numpy.array( cities3 )

    plt.plot(cities[:,0],cities[:,1],'bo-')
    #plt.scatter(cities[:,0], cities[:,1],s=50,c='k')

    if gain < 0:
        plt.scatter(cities2[:,0], cities2[:,1],c='r',s=180)
        plt.plot(cities2[:,0],cities2[:,1],c='r',linewidth=2)
        plt.scatter(cities3[:,0], cities3[:,1],c='b',s=150)
        plt.plot(cities3[:,0],cities3[:,1],c='r',linewidth=2)

    else:
        plt.scatter(cities2[:,0], cities2[:,1],c='g',s=180)
        plt.plot(cities2[:,0],cities2[:,1],c='g',linewidth=2)
        plt.scatter(cities3[:,0], cities3[:,1],c='b',s=150)
        plt.plot(cities3[:,0],cities3[:,1],c='g',linewidth=2)

    plt.axis( [-100,4100,-100,2100] )
    plt.axis('off')

    plt.title('(3)  2-Opt Tour {:6.1f}'.format(l_min))
    plt.savefig( ("%05d" % pic)+'.png')
    plt.clf()
    pic += 1
    print pic



#-----------------------------------------------------------------------------
#
#    Before 2opt             After 2opt
#       Y   Z                    Y   Z
#       O   O----->              O-->O---->
#      / \  ^                     \
#     /   \ |                      \
#    /     \|                       \
# ->O       O              ->O------>O
#   C       X                C       X
#
# In a 2opt optimization step we consider two nodes, Y and X.  (Between Y
# and X there might be many more nodes, but they don't matter.) We also
# consider the node C following Y and the node Z following X. i
#
# For the optimization we see replacing the edges CY and XZ with the edges CX
# and YZ reduces the length of the path  C -> Z.  For this we only need to
# look at |CY|, |XZ|, |CX| and |YZ|.   |YX| is the same in both
# configurations.
#
# If there is a length reduction we swap the edges AND reverse the direction
# of the edges between Y and X.
#
# In the following function we compute the amount of reduction in length
# (gain) for all combinations of nodes (X,Y) and do the swap for the
# combination that gave the best gain.
#

def optimize2opt(nodes, solution, number_of_nodes):
    best = 0
    best_move = None
    # For all combinations of the nodes
    for ci in range(0, number_of_nodes):
        for xi in range(0, number_of_nodes):
            yi = (ci + 1) % number_of_nodes  # C is the node before Y
            zi = (xi + 1) % number_of_nodes  # Z is the node after X

            c = solution[ ci ]
            y = solution[ yi ]
            x = solution[ xi ]
            z = solution[ zi ]
            # Compute the lengths of the four edges.
            cy = length( c, y )
            xz = length( x, z )
            cx = length( c, x )
            yz = length( y, z )

            # Only makes sense if all nodes are distinct
            if xi != ci and xi != yi:
                # What will be the reduction in length.
                gain = (cy + xz) - (cx + yz)
                # Is is any better then best one sofar?
                if gain > best:
                    # Yup, remember the nodes involved
                    best_move = (ci,yi,xi,zi)
                    best = gain

    print best_move, best
    if best_move is not None:
        (ci,yi,xi,zi) = best_move
        # This four are needed for the animation later on.
        c = solution[ ci ]
        y = solution[ yi ]
        x = solution[ xi ]
        z = solution[ zi ]

        # Create an empty solution
        new_solution = range(0,number_of_nodes)
        # In the new solution C is the first node.
        # This we we only need two copy loops instead of three.
        new_solution[0] = solution[ci]

        n = 1
        # Copy all nodes between X and Y including X and Y
        # in reverse direction to the new solution
        while xi != yi:
            new_solution[n] = solution[xi]
            n = n + 1
            xi = (xi-1)%number_of_nodes
        new_solution[n] = solution[yi]

        n = n + 1
        # Copy all the nodes between Z and C in normal direction.
        while zi != ci:
            new_solution[n] = solution[zi]
            n = n + 1
            zi = (zi+1)%number_of_nodes
        # Create a new animation frame
        frame4(nodes, new_solution, number_of_nodes, c, y, x, z, gain)
        return (True,new_solution)
    else:
        return (False,solution)


#-----------------------------------------------------------------------------
# This is an SA optimization step.
# It uses the same principle as the optimize2opt with the following
# differences:
#
# (1) Instead of all combinations of (X,Y) is picks a single combination
# at random.
#
# (1) Instead of only doing an edge swap if it reduces the length, it
# sometimes (depending on chance) also does a swap that INCREASES the length.
# How often this happens depends on the temperature t and the gain.
# For high temperatures this happens often and large negative gains are accepted,
# but the lower the temperature the less often it happens and only small
# negative gains are accepted.
#

def sa_optimize_step(nodes, solution, number_of_nodes, t):
    global nn
    # Pick X and Y at random.
    ci = random.randint(0, number_of_nodes-1)
    yi = (ci + 1) % number_of_nodes
    xi = random.randint(0, number_of_nodes-1)
    zi = (xi + 1) % number_of_nodes

    if xi != ci and xi != yi:
        c = solution[ci]
        y = solution[yi]
        x = solution[xi]
        z = solution[zi]
        cy = length(c, y)
        xz = length(x, z)
        cx = length(c, x)
        yz = length(y, z)

        gain = (cy + xz) - (cx + yz)
        if gain < 0:
            # We only accept a negative gain conditionally
            # The probability is based on the magnitude of the gain
            # and the temperature.
            u = math.exp( gain / t )
        elif gain > 0.05:
            u = 1 # always except a good gain.
        else:
            u = 0 # No idea why I did this....

        # random chance, picks a number in [0,1)
        if (random.random() < u):
            nn = nn + 1
            #print "      ", gain
            # Make a new solution with both edges swapped.
            new_solution = range(0,number_of_nodes)
            new_solution[0] = solution[ci]
            n = 1
            while xi != yi:
                new_solution[n] = solution[xi]
                n = n + 1
                xi = (xi-1)%number_of_nodes
            new_solution[n] = solution[yi]
            n = n + 1
            while zi != ci:
                new_solution[n] = solution[zi]
                n = n + 1
                zi = (zi+1)%number_of_nodes

            # Create an animation frame for this step
            frame(nodes, new_solution, number_of_nodes, t, c, y, x, z, gain)

            return new_solution
        else:
            return solution
    else:
        return solution


#----------------------------------------------------------------------------
def greedy_algorithm(nodes):
    # Greedy Algorithm
    print 'Computing greedy path'

    free_nodes = nodes[:]
    solution = []
    n = free_nodes[0]
    free_nodes.remove(n)
    solution.append( n )
    while len(free_nodes) > 0:
        print(len(free_nodes))
        min_l = None
        min_n = None
        for c in free_nodes:
            l = length( c, n )
            if min_l is None:
                min_l = l
                min_n = c
            elif l < min_l:
                min_l = l
                min_n = c
        solution.append(min_n)
        free_nodes.remove(min_n)
        n = min_n
    return solution

#-----------------------------------------------------------------------------
def two_opt_algorithm(nodes, number_of_nodes):
    # Create an initial solution
    solution = [n for n in nodes]
    go = True
    # Try to optimize the solution with 2opt until
    # no further optimization is possible.
    while go:
        (go,solution) = optimize2opt(nodes, solution, number_of_nodes)
    return solution

#-----------------------------------------------------------------------------
def sa_algorithm(nodes, number_of_nodes):
    # Create an initial solution that we can improve upon.
    solution = [n for n in nodes]

    # The temperature t. This is the most important parameter of the SA
    # algorithm. It starts at a high temperature and is then slowly decreased.
    # Both rate of decrease and initial values are parameters that need to be
    # tuned to get a good solution.

    # The initial temperature.  This should be high enough to allow the
    # algorithm to explore many sections of the search space.  Set too high it
    # will waste a lot of computation time randomly bouncing around the search
    # space.
    t = 100

    # Length of the best solution so far.
    l_min = total_length( nodes, solution )
    best_solution = []
    i = 0
    while t > 0.1:
        i = i + 1
        # Given a solution we create a new solution
        solution = sa_optimize_step(nodes, solution, number_of_nodes, t)
        # every ~200 steps
        if i >= 200:
            i = 0
            # Compute the length of the solution
            l = total_length( nodes, solution )
            print "    ", l, t, nn
            # Lower the temperature.
            # The slower we do this, the better then final solution
            # but also the more times it takes.
            t = t*0.9995

            # See if current solution is a better solution then the previous
            # best one.
            if l_min is None: # TODO: This can be removed, as l_min is set above.
                l_min = l
            elif l < l_min:
                # Yup it is, remember it.
                l_min = l
                print "++", l, t
                best_solution = solution[:]
            else:
                pass

    return best_solution


# From: http://stackoverflow.com/questions/16625507/python-checking-if-point-is-inside-a-polygon
# (Patrick Jordan)
# Modified to work with nodes.
def point_in_poly(x, y, poly):
    n = len(poly)
    inside = False

    p1x = poly[0].x
    p1y = poly[0].y
    for i in range(n+1):
        p2x = poly[i % n].x
        p2y = poly[i % n].y
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def area_triangle(n1, n2, n3):
    # Area of triangle via Heron's Formula
    # <https://en.wikipedia.org/wiki/Heron%27s_formula>
    a = length(n1, n2)
    b = length(n2, n3)
    c = length(n3, n1)
    p = (a + b + c)/2.0
    area = math.sqrt(p*(p-a)*(p-b)*(p-c))
    return area

def miss_perry_s_compass(nodes, number_of_nodes):
    # Compute Center of all nodes
    sum_x = 0
    sum_y = 0
    for n in nodes:
        sum_x += n.x
        sum_y += n.y
    c_x = sum_x / number_of_nodes
    c_y = sum_y / number_of_nodes
    # Make a node for this center.
    center_node = Node(-1, c_x, c_y)
    sorted_nodes = []
    done = [False] * number_of_nodes

    # Sort the nodes based on the distance from the center node.
    for i in range(number_of_nodes):
        max_l = -1
        furthest = None
        for j in range(number_of_nodes):
            if done[j]:
                pass
            else:
                l = length(center_node, nodes[j])
                if l > max_l:
                    furthest = j
                    max_l = l
        sorted_nodes.append(nodes[furthest])
        done[furthest] = True

    # Create initial polygon
    solution = [sorted_nodes[0], sorted_nodes[1], sorted_nodes[2]]
    for i in range(3, number_of_nodes):
        new_node = sorted_nodes[i]
        closest = sorted_nodes[0]
        min_l = length(closest, new_node)
        index_in_list = 0
        for j in range(1, i):
            l = length(sorted_nodes[j], new_node)
            if l < min_l:
                index_in_list = j
                closest = sorted_nodes[j]
                min_l = l

        # Is the node inside or outside the polygon?
        if point_in_poly(new_node.x, new_node.y, solution):
            idx_before = (index_in_list - 1) % i
            idx_after  = (index_in_list + 1) % i
            # it is Inside
            area1 = area_triangle(new_node, closest, solution[idx_before])
            area2 = area_triangle(new_node, closest, solution[idx_after])
            if area1 < area2:
                # Insert new node between closest and next
                pass
            else:
                # Insert 
                pass
            pass
        else:
            # it is outside
            pass

    return sorted_nodes

#-----------------------------------------------------------------------------
if False:
    # Experiment with Perry Algorithm
    do_greedy = False
    do_intro  = False
    do_perry  = True
    do_2opt   = False
    do_sa     = False
else:
    # Production
    do_greedy = True
    do_intro  = True
    do_perry  = False
    do_2opt   = True
    do_sa     = True

def create_animation(nodes):
    global nn
    global l_min
    number_of_nodes = len( nodes )
    print('Size {}'.format( number_of_nodes ))

    if do_greedy:
        # Greedy Algorithm
        print 'Computing greedy path'
        solution = greedy_algorithm(nodes)
    else:
        # For debugging
        solution = [n for n in nodes] 

    if do_intro:
        # Only cities
        solution0 = [n for n in nodes]
        for i in range(2, number_of_nodes):
            s = solution0[0:i]
            framec(s, number_of_nodes)
        # Show all cities for an additional 20 frames.
        for i in range(20):
            framec(s, number_of_nodes)

        # Animate the Random Search
        for i in range(2, number_of_nodes):
            s = solution0[0:i]
            frame0(s, nodes, total_length(nodes, s), "(1)  Random Path")
        s = solution0
        for i in range(60):
            frame0(s, nodes, total_length(nodes, s), "(1)  Random Path")

        # Animate the Greedy Search
        for i in range(2, number_of_nodes):
            s = solution[0:i]
            frame0(s, nodes, total_length(nodes, s), "(2)  Greedy Search")
        s = solution
        for i in range(60):
            frame0(s, nodes, total_length(nodes, s), "(2)  Greedy Search")

    # Under construction
    if do_perry:
        solution = miss_perry_s_compass(nodes, number_of_nodes)
        for i in range(2, number_of_nodes):
            s = solution[0:i]
            frame0(s, nodes, total_length(nodes, s), "(1)  Random Path")
        for i in range(60):
            frame0(solution, nodes, total_length(nodes, s), "(3)  Miss Perry")

    if do_2opt:
        print("2-Opt")
        # Run 2-Opt algorithm and create animation frames for each step
        s = two_opt_algorithm(nodes, number_of_nodes)
        # Show the best solution for an additional 60 frames.
        for i in range(60):
            frame0(s, nodes, total_length(nodes, s), "(4)  2-Opt")

    if do_sa:
        #=== Simulated Annealing 
        print("SA")
        # Run SA algorithm and create animation frames for each step
        s = sa_algorithm(nodes, number_of_nodes)
        # Show the best solution for an additional 60 frames.
        for i in range(60):
            frame0(s, nodes, total_length(nodes, s), "(5)  SA")

    return s

#-----------------------------------------------------------------------------
def read_problem(problem_file_name):
    nodes = []
    with open(problem_file_name) as inpf:
        first_line = inpf.readline()
        node_count = int(first_line)
        i = 0
        for line in inpf:
            parts = line.split()
            nodes.append(Node(i, float(parts[0]), float(parts[1])))
            i = i + 1

    return nodes

#-----------------------------------------------------------------------------
def solve(problem_file_name):
    # This it to make sure we get the same answer each time.
    random.seed(8111142)
    solution_string = None
    nodes = read_problem(problem_file_name)

    solution = create_animation(nodes)

    objective = total_length(nodes, solution)
    solution_string = str(objective) + ' 0\n'
    solution_string += ' '.join(map(lambda x: str(x.id), solution))

    return solution_string


if __name__ == '__main__':
    print solve('problem3.dat')

