import CoolProp.CoolProp as CP #I love coolprop! ~ units: http://www.coolprop.org/v4/apidoc/CoolProp.html


v  = 1/CP.PropsSI('D', 'T', 275, 'P', 52e5, "Nitrogen")
m = 0.12

V = v*m
print(V)