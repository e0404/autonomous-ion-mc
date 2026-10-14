import json,sys
t=open(sys.argv[1]).read()
j=t.split("#JSON-BEGIN")[1].split("#JSON-END")[0] if "#JSON-BEGIN" in t else t
d=json.loads(j)
def walk(o,p=""):
    if isinstance(o,dict):
        for k,v in o.items(): walk(v,p+"/"+k)
    elif isinstance(o,list):
        if len(o)<12 and all(not isinstance(x,(dict,list)) for x in o): print(p,o)
        else:
            for i,x in enumerate(o[:int(sys.argv[2]) if len(sys.argv)>2 else 6]): walk(x,p+f"[{i}]")
    else: print(p,o)
walk(d)
