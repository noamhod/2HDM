import matrix2py
import argparse
parser = argparse.ArgumentParser(description='run matrix !')
#########
parser.add_argument('-path',  metavar='--path',  required=True, help='library path')
#########
parser.add_argument('-g1e',  metavar='gluon[1].e',  required=True, help='g[1].e')
parser.add_argument('-g1px', metavar='gluon[1].px', required=True, help='g[1].px')
parser.add_argument('-g1py', metavar='gluon[1].py', required=True, help='g[1].py')
parser.add_argument('-g1pz', metavar='gluon[1].pz', required=True, help='g[1].pz')
#########                             
parser.add_argument('-g2e',  metavar='gluon[2].e',  required=True, help='g[2].e')
parser.add_argument('-g2px', metavar='gluon[2].px', required=True, help='g[2].px')
parser.add_argument('-g2py', metavar='gluon[2].py', required=True, help='g[2].py')
parser.add_argument('-g2pz', metavar='gluon[2].pz', required=True, help='g[2].pz')
#########                             
parser.add_argument('-t1e',  metavar='top[1].e',  required=True, help='t[1].e')
parser.add_argument('-t1px', metavar='top[1].px', required=True, help='t[1].px')
parser.add_argument('-t1py', metavar='top[1].py', required=True, help='t[1].py')
parser.add_argument('-t1pz', metavar='top[1].pz', required=True, help='t[1].pz')
#########                             
parser.add_argument('-t2e',  metavar='top[2].e',  required=True, help='t[2].e')
parser.add_argument('-t2px', metavar='top[2].px', required=True, help='t[2].px')
parser.add_argument('-t2py', metavar='top[2].py', required=True, help='t[2].py')
parser.add_argument('-t2pz', metavar='top[2].pz', required=True, help='t[2].pz')
#########
args = parser.parse_args()

def invert_momenta(p):
   """ fortran/C-python do not order table in the same order"""
   new_p = []
   for i in range(len(p[0])):  new_p.append([0]*len(p))
   for i, onep in enumerate(p):
      for j, x in enumerate(onep):
         new_p[j][i] = x
   return new_p

matrix2py.initialise(args.path+'/param_card.dat')
p = [[ args.g1e, args.g1px, args.g1py, args.g1pz ],
     [ args.g2e, args.g2px, args.g2py, args.g2pz ],
     [ args.t1e, args.t1px, args.t1py, args.t1pz ],
     [ args.t2e, args.t2px, args.t2py, args.t2pz ]]
P=invert_momenta(p)

### running with:
# python matrix.py -g1e 0.5000000E+03 -g1px 0.0000000E+00 -g1py 0.0000000E+00 -g1pz 0.5000000E+03 -g2e 0.5000000E+03 -g2px 0.0000000E+00 -g2py 0.0000000E+00 -g2pz -500 -t1e 0.5000000E+03 -t1px 0.1109243E+03 -t1py 0.4448308E+03 -t1pz -199.5529 -t2e 0.5000000E+03 -t2px -110.9243 -t2py -444.8308 -t2pz  0.1995529E+03
# p = [[   0.5000000E+03,  0.0000000E+00,  0.0000000E+00,  0.5000000E+03],
#     [   0.5000000E+03,  0.0000000E+00,  0.0000000E+00, -0.5000000E+03],
#     [   0.5000000E+03,  0.1109243E+03,  0.4448308E+03, -0.1995529E+03],
#     [   0.5000000E+03, -0.1109243E+03, -0.4448308E+03,  0.1995529E+03]]

alphas = 0.13
nhel   = 0 # means sum over all helicity                                                                                                                  
me2 = matrix2py.get_me(P, alphas, nhel)
print me2
