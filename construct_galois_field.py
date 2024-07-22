from zp_operations import GaloisField

p = 3
f = [1, 0, 1]

gf = GaloisField(p, f)
gf.print_table(operation="addition", mode="polynomial")
gf.print_table(operation="addition", mode="integer")

gf.print_table(operation="multiplication", mode="polynomial")
gf.print_table(operation="multiplication", mode="integer")
