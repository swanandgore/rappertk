import os

def main() :
    print "delete all"
    #print "set stick_radius, 0.5"
    print "bg_color grey"

    print "load 121p_c.pdb"
    print "hide all"

    print "load 121p_1.map, format = ccp4"
    print "isosurface surf, 121p_1.map, 1"
    print "color red, surf"

    print "show cartoon"
    print "color green, 121p_c"
    #print "show spheres, name CA"
    #print "set sphere_transparency, 0.9"
    #print "set cartoon_transparency, 0.6"
    colors = ["red", "green", "blue", "magenta"]
    for i in range(0,1000) :
        clr = colors[ i%len(colors) ]
        fn = "snap%d.pdb" % i
        if not os.access(fn, os.R_OK) : continue
        print "load %s, snap" % fn
        print "hide everything, snap"
        #print "show ribbon, snap"
        print "show sticks, snap and (name N or name CA or name C or name O)"
        print "set all_states, 1"
        print "load 121p_c.pdb, x"
        print "zoom x and name CA"
        print "delete x"
        print  """set_view (\
                -0.013586010,   -0.327523857,   -0.944741488,\
                -0.998541236,   -0.044929720,    0.029940346,\
                -0.052251201,    0.943772912,   -0.326435059,\
                -0.000001247,    0.000001982, -147.574310303,\
                 0.137311488,   -1.178305864,   -0.623766601,\
               113.509162903,  181.639450073,    0.000000000 )"""
        print "ray"
        print "png b%d.png" % i
        print "delete snap"
    print "convert -adjoin b*png x.gif"

if __name__ == "__main__" : main()
