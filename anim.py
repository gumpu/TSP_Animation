#!/usr/bin/python

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
v = [-100,4100,-100,2100]

def framec( solution, nc ):
    global pic
    cities = [ (n.x,n.y) for n in solution ]
    cities = numpy.array( cities )
    plt.axis( [-100,4100,-100,2100] )
    plt.axis('off')
    plt.plot(cities[:,0],cities[:,1],'ko')
    plt.title('{} Cities, 4 Search Algorithms'.format(nc))
    plt.savefig( ("%05d" % pic)+'.png')
    plt.clf()
    print pic
    pic += 1

def frame0( solution, all, l, title ):
    global pic
    cities = [ (n.x,n.y) for n in solution ]
    cities.append( ( solution[0].x, solution[0].y ) )
    cities = numpy.array( cities )
    all = [ (n.x,n.y) for n in all ]
    all = numpy.array( all )

    plt.axis( [-100,4100,-100,2100] )
    plt.axis('off')
    plt.plot(cities[:,0],cities[:,1],'bo-')
    plt.plot(all[:,0],all[:,1],'ko')
    plt.title('{}  Tour length {:.1f}'.format(title,l))
    plt.savefig( ("%05d" % pic)+'.png')
    plt.clf()
    print pic
    pic += 1

nn = 0
def frame( nodes, solution, sn, t, c,y,x,z, gain ):
    global pic
    global nn

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

    plt.title('(4)  SA Temp {:4.1f} Best Tour {:6.1f}\nSwaps {}  Gain {:12.2f} '.format(t,l_min,nn,gain))
    plt.savefig( ("%05d" % pic)+'.png')
    plt.clf()
    pic += 1
    print pic

def frame4( nodes, solution, sn, t, c,y,x,z, gain ):
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
def create_animation( nodes ):
    global nn
    global l_min
    sn = len( nodes )
    print 'Size {}'.format( sn )

    if True:
        # Greedy Algorithm
        print 'Computing greedy path'

        free_nodes = nodes[:]
        solution = []
        n = free_nodes[0]
        free_nodes.remove(n)
        solution.append( n )
        while len( free_nodes ) > 0:
            print len( free_nodes ),
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
    else:
        solution = [ n for n in nodes ]

    if True:
        # Only cities
        solution0 = [ n for n in nodes ]
        for i in range(2,sn):
            s = solution0[0:i]
            framec(s,sn)
        for i in range(20):
            framec(s,sn)

        # Random Search
        for i in range(2,sn):
            s = solution0[0:i]
            frame0(s,solution0, total_length(nodes,s), "(1)  Random Path")
        s = solution0
        for i in range(60):
            frame0(s,solution, total_length(nodes,s), "(1)  Random Path")

        # Greedy
        for i in range(2,sn):
            s = solution[0:i]
            frame0(s,solution, total_length(nodes,s), "(2)  Greedy Search")
        s = solution
        for i in range(60):
            frame0(s,solution, total_length(nodes,s), "(2)  Greedy Search")


    print "2-Opt"
    solution = [ n for n in nodes ]
    t = 100
    go = True
    while go:
        (go,solution) = optimize2opt( nodes, solution, sn, t )
    s = solution
    for i in range(60):
        frame0(s,solution, total_length(nodes,s), "(3)  2-Opt")



    print "SA"
    solution = [ n for n in nodes ]
    t = 100
    l_min = total_length( nodes, solution )
    best_solution = []
    i = 0
    while True:
        i = i + 1
        solution = optimize( nodes, solution, sn, t )
        if i >= 200:
            i = 0
            l = total_length( nodes, solution )
            print "    ", l, t, nn
            t = t*0.9995

            if l_min is None:
                l_min = l
            elif l < l_min:
                l_min = l
                print "++", l, t
                best_solution = solution[:]
            else:
                pass
        if t < 0.1:
            break

    s = best_solution
    for i in range(60):
        frame0(s, solution, total_length(nodes,s), "(4)  SA")


    return best_solution

def optimize2opt( nodes, solution, sn, t ):
    best = 0
    best_move = None
    for ci in range(0, sn):
        for xi in range(0, sn):
            yi = (ci + 1) % sn
            zi = (xi + 1) % sn

            c = solution[ ci ]
            y = solution[ yi ]
            x = solution[ xi ]
            z = solution[ zi ]
            cy = length( c, y )
            xz = length( x, z )
            cx = length( c, x )
            yz = length( y, z )


            if xi != ci and xi != yi:
                gain = (cy + xz) - (cx + yz)
                if gain > best:
                    best_move = (ci,yi,xi,zi)
                    best = gain

    print best_move, best
    if best_move is not None:
        (ci,yi,xi,zi) = best_move
        c = solution[ ci ]
        y = solution[ yi ]
        x = solution[ xi ]
        z = solution[ zi ]

        new_solution = range(0,sn)
        new_solution[0] = solution[ci]
        n = 1
        while xi != yi:
            new_solution[n] = solution[xi]
            n = n + 1
            xi = (xi-1)%sn
        new_solution[n] = solution[yi]
        n = n + 1
        while zi != ci:
            new_solution[n] = solution[zi]
            n = n + 1
            zi = (zi+1)%sn
        frame4( nodes, new_solution, sn, t, c,y,x,z, gain )
        return (True,new_solution)
    else:
        return (False,solution)



def optimize( nodes, solution, sn, t ):
    global nn
    ci = random.randint(0, sn-1)
    yi = (ci + 1) % sn
    xi = random.randint(0, sn-1)
    zi = (xi + 1) % sn

    if xi != ci and xi != yi:
        c = solution[ ci ]
        y = solution[ yi ]
        x = solution[ xi ]
        z = solution[ zi ]
        cy = length( c, y )
        xz = length( x, z )
        cx = length( c, x )
        yz = length( y, z )

        gain = (cy + xz) - (cx + yz)
        if gain < 0:
            u = math.exp( gain / t )
        elif gain > 0.05:
            u = 1
        else:
            u = 0

        if (random.random() < u):
            nn = nn + 1
            #print "      ", gain
            new_solution = range(0,sn)
            new_solution[0] = solution[ci]
            n = 1
            while xi != yi:
                new_solution[n] = solution[xi]
                n = n + 1
                xi = (xi-1)%sn
            new_solution[n] = solution[yi]
            n = n + 1
            while zi != ci:
                new_solution[n] = solution[zi]
                n = n + 1
                zi = (zi+1)%sn

            frame( nodes, new_solution, sn, t, c,y,x,z, gain )

            return new_solution
        else:
            return solution
    else:
        return solution


def solve( input_data ):
    # This it to make sure we get the same anwser each time.
    random.seed(8111142)

    lines = input_data.split('\n')
    node_count = int(lines[0])
    nodes = []

    for i in range(1, node_count+1):
        line = lines[i]
        parts = line.split()
        nodes.append( Node( i-1, float(parts[0]), float(parts[1]) ) )

    solution = create_animation( nodes )

    objective = total_length( nodes, solution )
    solution_string = str( objective ) + ' 0\n'
    solution_string += ' '.join( map( lambda x: str(x.id), solution ) )
    return solution_string


if __name__ == '__main__':
    fileLocation = 'problem3.dat'
    input_data_file = open(fileLocation, 'r')
    input_data = ''.join(input_data_file.readlines())
    input_data_file.close()
    print solve( input_data )

