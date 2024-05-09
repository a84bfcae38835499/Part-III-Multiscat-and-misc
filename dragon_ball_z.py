from sympy.parsing.mathematica import parse_mathematica
from sympy import var

S,Theta,C,ratio = var('S,Theta,C,ratio')
expr='((1/(2 (-1+2 Theta)))(4 Theta-(16 Theta^2+(1/(Theta^2-2 Theta^3+Theta^4))^0.5 4 (-1+2 Theta) (-Theta+Theta^2-Theta^3+3 Theta^4-2 Theta^5-(Theta^2 ratio-2 Theta^3 ratio+Theta^4 ratio))^0.5)))^2'
my_expr_in_python  = parse_mathematica(expr)

S,Theta,C,ratio = var('S,Theta,C,ratio')
expr='-(C/(2 ratio^0.5 Theta))'
my_expr_in_python  = parse_mathematica(expr)
