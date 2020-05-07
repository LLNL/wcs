import pandas
 
stoic = pandas.read_csv('StoicMat_PARCDL2.csv',index_col=0)

import re

antimonyFile = open('BigModel_PARCDL2.txt','r')
constants = {}
for line in antimonyFile:
    m = re.search(r'([A-Za-z0-9_]+)\s*=\s*([-.0-9e+]+);', line)
    if m:
        constants[m.group(1)] = m.group(2)

antimonyFile = open('BigModel_PARCDL2.txt','r')
rxnRates = {}
for line in antimonyFile:
    m = re.search(r'(\w*)\s*:.*; ([^ ].*?)\s*;', line)
    if m:
        rxnRates[m.group(1)] = m.group(2)

dependentConstants = {}
dependentSpecies = {}
for (rxn,rate) in rxnRates.items():
    allVars = re.findall(r'([A-Za-z_][A-Za-z0-9_]*)',rate)
    thisDepSpecies = []
    thisDepConstants = []
    for var in allVars:
        if re.match(r'k.*|Cytoplasm|Extracellular|Nucleus|Mitochondrion',var):
            thisDepConstants.append(var)
        else:
            thisDepSpecies.append(var)
    dependentConstants[rxn] = thisDepConstants
    dependentSpecies[rxn] = thisDepSpecies

print("""
<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns">
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
        http://graphml.graphdrawing.org/xmlns/1.1/graphml.xsd">

    <!-- vertex (species and reactions) attributes -->
    <key id="v_label" for="node" attr.name="v_label" attr.type="string"/>
    <key id="v_type" for="node" attr.name="v_type" attr.type="int"/>
    <key id="s_count" for="node" attr.name="s_count" attr.type="int">
      <default>0</default>
    </key>
    <key id="r_const" for="node" attr.name="r_const" attr.type="double">
      <default>1.0</default>
    </key>
    <key id="r_rate" for="node" attr.name="r_rate" attr.type="string"/>

    <!-- edge attributes -->
    <key id="e_label" for="edge" attr.name="e_label" attr.type="string"/>
    <key id="e_stoic" for="edge" attr.name="e_stoic" attr.type="int">
        <default>1</default>
    </key>

    <graph id="simple" edgedefault="directed">
""")
print("""
        <!-- species -->""")

defaultVolume = 6.0221409e+23*float(constants['Cytoplasm'])/1e6;
if defaultVolume > 10000000:
    defaultVolume = defaultVolume / 1000

for speciesName in stoic.index.values:
    init_cnt = int(float(constants[speciesName.strip()])*defaultVolume)
    if init_cnt > 1000000000:
        init_cnt = init_cnt / 1000
    print("""
        <node id="s_%(name)s">
            <data key="v_label">%(name)s</data>
            <data key="v_type">1</data>
            <data key="s_count">%(init)d</data>
        </node>""" % {'name': speciesName.strip(), 'init': init_cnt})
print("""
        <!-- reactions -->""")

for (rxnName, rate) in rxnRates.items():
    actualRate = '\n'
    for ccc in set(dependentConstants[rxnName]):
        actualRate += 'var '+ccc+' := '+constants[ccc]+';\n'
    actualRate += 'm_rate := '+rate+';'
    print("""
        <node id="r_%(rxnName)s">
            <data key="v_label">%(rxnName)s</data>
            <data key="v_type">2</data>
            <data key="r_rate">%(rate)s</data>
        </node>""" % {"rxnName": rxnName, 'rate': actualRate})
if 1:
    for (rxnName, rxn) in stoic.iteritems():
        reactants = {}
        products = {}
        for (Cname,Cval) in rxn.items():
            Cclean = str(Cname).strip()
            if Cval <= -1:
                reactants[Cclean] = -Cval
            elif Cval >= 1:
                products[Cclean] = Cval
        count = 1
        for (reactant, sss) in reactants.items():
            print('''
            <edge id="r_%(rxn)s_r%(c)s" source="s_%(species)s" target="r_%(rxn)s">
               <data key="e_stoic">%(stoic)s</data>
            </edge>''' % {'rxn':rxnName, 'species': reactant, 'c': count, 'stoic': sss})
            count += 1
        count = 1
        for (product, sss) in products.items():
            print('''
            <edge id="r_%(rxn)s_p%(c)s" source="r_%(rxn)s" target="s_%(species)s">
               <data key="e_stoic">%(stoic)s</data>
            </edge>''' % {'rxn':rxnName, 'species': product, 'c': count, 'stoic': sss})
            count += 1
        count = 1
        ddd = set(dependentSpecies[rxnName]) - set(reactants.keys())
        for dep in ddd:
            print('''
            <edge id="r_%(rxn)s_d%(c)s" source="s_%(species)s" target="r_%(rxn)s">
               <data key="e_stoic">0</data>
            </edge>''' % {'rxn':rxnName, 'species': dep, 'c': count})
            count += 1

print("""

    </graph>
</graphml>""")
