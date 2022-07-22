import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import itertools


def pol2cart(r, theta):
    z = r * np.exp(1j * theta)
    x, y = z.real, z.imag
    return x, y


def cart2pol(x, y):
    z = x + y * 1j
    r, theta = np.abs(z), np.angle(z)
    return r, theta


class Node:
    nodeNumber = 0
    names = set()

    @property
    def name(self):
        return self._name

    def __init__(self, type: int, voltage: float, theta: float, PGi: float, QGi: float, PLi: float, QLi: float,
                 Qmin: float, Qmax: float, name=None):
        self.nodeNumber = Node.nodeNumber
        Node.nodeNumber += 1
        self.voltage = voltage
        self.type = type
        self.theta = theta
        self.PGi = PGi
        self.QGi = QGi
        self.PLi = PLi
        self.QLi = QLi
        self.Qmin = Qmin
        self.Qmax = Qmax
        self.vLf = self.voltage
        self.thetaLf = self.theta

        if name in Node.names:
            Node.nodeNumber -= 1
            raise NameError("Already used name '%s'." % name)
        if name is None:
            self._name = str(self.nodeNumber)
        else:
            self._name = name

        Node.names.add(self.name)

    @property
    def vm(self):
        return self.voltage * np.exp(self.theta * 1j)

    @property
    def vmLf(self):
        return self.vLf * np.exp(self.thetaLf * 1j)


class Line:
    lineNumber = 0
    names = set()

    @property
    def name(self):
        return self._name

    def __init__(self, fromNode: Node, toNode: Node, r: float, x: float, b_half: float, x_prime: float, name=None):
        Line.lineNumber += 1
        self.lineNumber = Line.lineNumber
        self.fromNode = fromNode
        self.toNode = toNode
        self.r = r
        self.x = x
        self.b_half = b_half
        self.x_prime = x_prime
        self.z = self.r + self.x * 1j
        self.y = 1 / self.z
        self.b = self.b_half * 1j

        # avoiding the same name for two different lines
        if name in Node.names:
            Line.lineNumber -= 1
            raise NameError("Already used name '%s'." % name)

        if name is None:
            self._name = str(self.lineNumber)
        else:
            self._name = name
        Line.names.add(self.name)


class Grid:
    def __init__(self, nodes: list, lines: list):
        self.nodes = nodes
        self.lines = lines
        self.Y = np.zeros((self.nb, self.nb), dtype=complex)
        self.nl = len(self.lines)
        self.create_matrix()
        self.Pl = np.vstack([node.PLi for node in self.nodes])
        self.Ql = np.vstack([node.QLi for node in self.nodes])
        self.Pg = np.vstack([node.PGi for node in self.nodes])
        self.Qg = np.vstack([node.QGi for node in self.nodes])
        self.Psp = self.Pg - self.Pl
        self.Qsp = self.Qg - self.Ql

    @property
    def nb(self):
        fromBus = [line.fromNode.nodeNumber for line in self.lines]
        toBus = [line.toNode.nodeNumber for line in self.lines]
        return max(max(fromBus), max(toBus)) + 1 # +1 because of the 0-indexing

    def get_node_by_number(self, number: int):
        for node in self.nodes:
            if node.nodeNumber == number:
                return node
        raise NameError("No node with number %d." % number)

    def get_line_by_number(self, number: int):
        for line in self.lines:
            if line.lineNumber == number:
                return line
        raise NameError("No line with number %d." % number)

    def get_lines_by_node(self, nodeNumber):
        lines = [line for line in self.lines if
                 (line.toNode.nodeNumber == nodeNumber or line.fromNode.nodeNumber == nodeNumber)]
        return lines

    @property
    def pq_nodes(self):
        pq_nodes = [node for node in self.nodes if node.type == 3]
        return pq_nodes

    @property
    def pv_nodes(self):
        pq_nodes = [node for node in self.nodes if node.type == 2]
        return pq_nodes

    def create_matrix(self):
        # off diagonal elements
        for k in range(self.nl):
            line = self.lines[k]
            fromNode = line.fromNode.nodeNumber
            toNode = line.toNode.nodeNumber
            self.Y[fromNode, toNode] -= line.y/line.x_prime
            self.Y[toNode, fromNode] = self.Y[fromNode, toNode]

        # diagonal elements
        for m in range(self.nb):
            for n in range(self.nl):
                line = self.lines[n]
                if line.fromNode.nodeNumber == m:
                    self.Y[m, m] += line.y/(line.x_prime**2) + line.b
                elif line.toNode.nodeNumber == m:
                    self.Y[m, m] += line.y+ line.b

    def calculateLf(self, BMva=100):
        Vm = np.vstack([node.vmLf for node in self.nodes]).reshape((self.nb, -1))
        self.I = np.matmul(self.Y, Vm)
        Iij = np.zeros((self.nb, self.nb), dtype=complex)
        Sij = np.zeros((self.nb, self.nb), dtype=complex)

        self.Im = abs(self.I)
        self.Ia = np.angle(self.I)

        for node in self.nodes:
            m = node.nodeNumber  # node index
            lines = self.get_lines_by_node(nodeNumber=node.nodeNumber)
            for line in lines:
                if line.fromNode.nodeNumber == m:
                    p = line.toNode.nodeNumber  # index to
                    if m != p:
                        Iij[m,p] = -(line.fromNode.vmLf - line.toNode.vmLf * line.x_prime) * self.Y[m,p]/(line.x_prime ** 2) + line.b_half / (line.x_prime ** 2) * line.fromNode.vmLf
                        Iij[p,m] = - (line.toNode.vmLf - line.fromNode.vmLf / line.x_prime) * self.Y[p,m] + line.b_half * line.toNode.vmLf
                else:
                    p = line.fromNode.nodeNumber  # index from
                    if m != p:
                        Iij[m,p] = - (line.toNode.vmLf - line.fromNode.vmLf / line.x_prime) * self.Y[p,m] + line.b_half * line.toNode.vmLf
                        Iij[p,m] = - (line.fromNode.vmLf - line.toNode.vmLf) * self.Y[m,p] / (
                                    line.x_prime ** 2) + line.b_half / (line.x_prime ** 2) * line.fromNode.vmLf

        self.Iij = Iij
        self.Iijr = np.real(Iij)
        self.Iiji = np.imag(Iij)

        # line powerflows
        for m in range(self.nb):
            for n in range(self.nb):
                if n != m:
                    Sij[m,n] = self.nodes[m].vmLf * np.conj(self.Iij[m, n]) * BMva

        self.Sij = Sij
        self.Pij = np.real(Sij)
        self.Qij = np.imag(Sij)

        # line losses
        Lij = np.zeros(self.nl, dtype=complex)
        for line in self.lines:
            m = line.lineNumber - 1
            p = line.fromNode.nodeNumber
            q = line.toNode.nodeNumber
            Lij[m] = Sij[p, q] + Sij[q, p]

        self.Lij = Lij
        self.Lpij = np.real(Lij)
        self.Lqij = np.imag(Lij)

        # Bus power injection
        Si = np.zeros(self.nb, dtype=complex)
        for i in range(self.nb):
            for k in range(self.nb):
                Si[i] += np.conj(self.nodes[i].vmLf) * self.nodes[k].vmLf * self.Y[i, k] * BMva

        self.Si = Si
        self.Pi = np.real(Si)
        self.Qi = -np.imag(Si)
        self.Pg = self.Pi.reshape([-1,1]) + self.Pl.reshape([-1,1])
        self.Qg = self.Qi.reshape([-1,1]) + self.Ql.reshape([-1,1])

    def loadflow(self, tol=1, maxIter=10000, BMva=100):
        self.iter = 0
        Pg = self.Pg / BMva
        Qg = self.Qg / BMva
        Pl = self.Pl / BMva
        Ql = self.Ql / BMva
        Psp = self.Psp / BMva
        Qsp = self.Qsp / BMva
        G = np.real(self.Y)
        B = np.imag(self.Y)
        angles = np.zeros((self.nb,1))
        npv = len(self.pv_nodes)
        npq = len(self.pq_nodes)
        self.tolerances = []

        while self.iter < 20 or (tol > 1e-5 and self.iter < maxIter):
            self.iter += 1
            P = np.zeros((self.nb, 1))
            Q = np.zeros((self.nb, 1))

            #calculate P and Q
            for node in self.nodes:
                i = node.nodeNumber
                for k in range(self.nb):
                    P[i] += node.vLf*self.nodes[k].vLf*(G[i, k]*np.cos(angles[i]-angles[k]) + B[i,k]*np.sin(angles[i]-angles[k]))
                    Q[i] += node.vLf*self.nodes[k].vLf*(G[i, k]*np.sin(angles[i]-angles[k]) - B[i,k]*np.cos(angles[i]-angles[k]))
            self.P = P

            #calculate Q-limit violations
            if self.iter > 2 and self.iter <= 7:
                for n in range(1, self.nb):
                    if self.nodes[n].type == 2:
                        QG = Q[n] + Ql[n]
                        if QG < self.nodes[n].Qmin/BMva:
                            self.nodes[n].vLf += 0.01
                        elif QG > self.nodes[n].Qmax/BMva:
                            self.nodes[n].vLf -= 0.01

            #calculate changes in specified active and reactive power
            dPa = Psp - P
            dQa = Qsp - Q
            k = 0
            dQ = np.zeros((npq, 1))
            for node in self.pq_nodes:
                i = node.nodeNumber
                if node.type == 3:
                    dQ[k] = dQa[i]
                    k += 1
            dP = dPa[1:self.nb]
            M = np.vstack((dP, dQ))

            #calculate Jacobian. #
            #J1 is the derivative of P with respect to angles
            J1 = np.zeros((self.nb-1, self.nb-1))
            for i in range(self.nb-1):
                m = i + 1
                for k in range(self.nb-1):
                    n = k + 1
                    if n == m:
                        for n in range(self.nb):
                            J1[i, k] += self.nodes[m].vLf*self.nodes[n].vLf*(-G[m,n]*np.sin(angles[m]-angles[n]) + B[m,n]*np.cos(angles[m]-angles[n]))
                        J1[i,k] += -self.nodes[m].vLf**2*B[m,m]
                    else:
                        J1[i,k] = self.nodes[m].vLf*self.nodes[n].vLf*(G[m,n]*np.sin(angles[m]-angles[n]) - B[m,n]*np.cos(angles[m]-angles[n]))
            self.J1 = J1

            #J2 is the derivative of P with respect to V
            J2 = np.zeros((self.nb-1, npq))
            for i in range(self.nb-1):
                m = i + 1
                for k in range(npq):
                    n = self.pq_nodes[k].nodeNumber
                    if n == m:
                        for n in range(self.nb):
                            J2[i, k] += self.nodes[n].vLf*(G[m,n]*np.cos(angles[m]-angles[n]) + B[m,n]*np.sin(angles[m]-angles[n]))
                        J2[i,k] += self.nodes[m].vLf*G[m,m]
                    else:
                        J2[i,k] = self.nodes[m].vLf*(G[m,n]*np.cos(angles[m]-angles[n]) + B[m,n]*np.sin(angles[m]-angles[n]))
            self.J2 = J2


            #J3 is the derivative of Q with respect to angles
            J3 = np.zeros((npq, self.nb-1))
            for i in range(npq):
                m = self.pq_nodes[i].nodeNumber
                for k in range(self.nb-1):
                    n = k + 1
                    if n == m:
                        for n in range(self.nb):
                            J3[i, k] += self.nodes[m].vLf*self.nodes[n].vLf*(G[m,n]*np.cos(angles[m]-angles[n]) + B[m,n]*np.sin(angles[m]-angles[n]))
                        J3[i,k] += -self.nodes[m].vLf**2*G[m,m]
                    else:
                        J3[i,k] = self.nodes[m].vLf*self.nodes[n].vLf*(-G[m,n]*np.cos(angles[m]-angles[n]) - B[m,n]*np.sin(angles[m]-angles[n]))
            self.J3 = J3

            #J4 is the derivative of Q with respect to V
            J4 = np.zeros((npq, npq))
            for i in range(npq):
                m = self.pq_nodes[i].nodeNumber
                for k in range(npq):
                    n = self.pq_nodes[k].nodeNumber
                    if n == m:
                        for n in range(self.nb):
                            J4[i, k] += self.nodes[n].vLf*(G[m,n]*np.sin(angles[m]-angles[n]) - B[m,n]*np.cos(angles[m]-angles[n]))
                        J4[i,k] += -self.nodes[m].vLf*B[m,m]
                    else:
                        J4[i,k] = self.nodes[m].vLf*(G[m,n]*np.sin(angles[m]-angles[n]) - B[m,n]*np.cos(angles[m]-angles[n]))
            self.J4 = J4

            self.J = np.vstack((np.hstack((J1, J2)), np.hstack((J3, J4))))
            #end of Jacobian calculation
            #J X = M -> X = J^-1 M
            X = np.linalg.solve(self.J, M)
            dTh = X[0:self.nb-1]
            dV = X[self.nb-1:]

            #update Angles and Voltages
            angles[1:] += dTh #angles[0] is the angle of the slack bus
            k=0
            for i in range(1, self.nb):
                if self.nodes[i].type == 3:
                    self.nodes[i].vLf += dV[k].item()
                    k += 1
                self.nodes[i].thetaLf = angles[i].item()

            tol = max(abs(M))
            self.tolerances.append((tol))
            self.voltageLf = [self.nodes[i].vLf for i in range(self.nb)]
            self.thetaLf = [self.nodes[i].thetaLf for i in range(self.nb)]

        #the iteration is over; calculate the power flow
        self.calculateLf()

    def printResults(self):
        print("Newton Raphson Results:")
        print()
        print('| Bus |    V     |  Angle   |      Injection      |      Generation     |          Load      |')
        print('| No  |    pu    |  Degree  |     MW   |   MVar   |     MW   |  Mvar    |     MW  |     MVar |')
        for i in range(self.nb):
            print('| %3g | %8.3f | %8.3f | %8.3f | %8.3f | %8.3f | %8.3f |%8.3f | %8.3f |' % (i, self.nodes[i].vLf, self.nodes[i].thetaLf, self.Pi[i], self.Qi[i], self.Pg[i], self.Qg[i], self.Pl[i], self.Ql[i]))

        print('----------------------------------------------------------------------------------------------')
        print()
        print("Line flows and losses:")
        print()
        print('|From  |To    |     P    |     Q    | From | To   |    P     |    Q     |     Line Loss       |')
        print('|Bus   |Bus   |    MW    |    MVar  | Bus  | Bus  |    MW    |   MVar   |     MW   |    MVar  |')
        for i in range(self.nl):
            p = self.lines[i].fromNode.nodeNumber
            q = self.lines[i].toNode.nodeNumber
            print('| %4g | %4g | %8.2f | %8.2f | %4g | %4g | %8.2f | %8.2f | %8.2f | %8.2f |' % (p, q, self.Pij[p,q], self.Qij[p,q], q, p, self.Pij[q,p], self.Qij[q,p], self.Lpij[i], self.Lqij[i]))
        print('----------------------------------------------------------------------------------------------')
        print()
        print('Total active losses: {active_power:.2f}, Total reactive losses: {reactive_power:.2f}'.format(active_power=sum(self.Lpij).item(), reactive_power=sum(self.Lqij).item()))
