from Utilities.ElectricalObjects import *

def main():
    n1 = Node(1, 1.05, 0, 0.0, 0, 0, 0, 0, 0)
    n2 = Node(2, 1.05, 0, 0.5, 0, 0, 0, -0.5, 1.0)
    n3 = Node(2, 1.07, 0, 0.6, 0, 0, 0, -0.5, 1.5)
    n4 = Node(3, 1.0, 0, 0.0, 0, 0.7, 0.7, 0, 0)
    n5 = Node(3, 1.0, 0, 0.0, 0, 0.7, 0.7, 0, 0)
    n6 = Node(3, 1.0, 0, 0.0, 0, 0.7, 0.7, 0, 0)

    l1 = Line(n1, n2, 0.1, 0.2, 0.02, 1)
    l2 = Line(n1, n4, 0.05, 0.2, 0.02, 1)
    l3 = Line(n1, n5, 0.08, 0.3, 0.03, 1)
    l4 = Line(n2, n3, 0.05, 0.25, 0.03, 1)
    l5 = Line(n2, n4, 0.05, 0.1, 0.01, 1)
    l6 = Line(n2, n5, 0.1, 0.3, 0.02, 1)
    l7 = Line(n2, n6, 0.07, 0.2, 0.025, 1)
    l8 = Line(n3, n5, 0.12, 0.26, 0.025, 1)
    l9 = Line(n3, n6, 0.02, 0.1, 0.01, 1)
    l10 = Line(n4, n5, 0.2, 0.4, 0.04, 1)
    l11 = Line(n5, n6, 0.1, 0.3, 0.03, 1)

    nodes = [n1, n2, n3, n4, n5, n6]
    lines = [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11]

    grid = Grid(nodes=nodes, lines=lines)
    grid.loadflow()

if __name__ == '__main__':
    main()


