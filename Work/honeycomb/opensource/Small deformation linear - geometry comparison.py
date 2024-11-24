# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from odbAccess import *
import numpy as np
import string
import math

##Adding data to csv
topreactionvtxtlist=[]
bottomreactionvtxtlist=[]
topreactionrelativethicknesstxtlist=[]
bottomreactionrelativethicknesstxtlist=[]

#Lithitated Thickness (It has to be bigger than seedsize, or minimum seedsize ?) Doesn't work at 0.2

top_thickness= 0.4
bottom_thickness = 0.4
# Time step ( % of a second per increment, default a full second per increment)
timestep = 0.01


#Model number 
x= 7
eigenvalue= 1000

while eigenvalue > 1 :
    x+=1
    #LET'S CREATE MODEL#
    
    mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-%d' %(x))
    
    ##sketch geometry in nm
    #Bare in mind sheet size can be specified below 
    height = 250
    length = 11600
    
   ##Sketch (fixed)##
    mdb.models['Model-%d' %(x)].ConstrainedSketch(name='__profile__', sheetSize= 100000)
    mdb.models['Model-%d' %(x)].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
        point2=(length, height))
    mdb.models['Model-%d' %(x)].Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-%d' %(x)].parts['Part-1'].BaseShell(sketch=
        mdb.models['Model-%d' %(x)].sketches['__profile__'])
    del mdb.models['Model-%d' %(x)].sketches['__profile__']

    ##Partition Sketch (Main modify area)##

    #parameters for partition
    #Lithiated Silicon thickness Max must be < (height/2)
    bottomcornerx = 0
    bottomcornery = 0
    #Sketch setup and sketch origin transformation(if needed)
    mdb.models['Model-%d' %(x)].ConstrainedSketch(gridSpacing=1, name='__profile__', 
        sheetSize=1000, transform=
        mdb.models['Model-%d' %(x)].parts['Part-1'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-%d' %(x)].parts['Part-1'].faces.findAt((
        length-0.001,height-0.001,0),), sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0,0,0)))
    mdb.models['Model-%d' %(x)].parts['Part-1'].projectReferencesOntoSketch(filter=
        COPLANAR_EDGES, sketch=mdb.models['Model-%d' %(x)].sketches['__profile__'])
    #Partition top sectoin
    mdb.models['Model-%d' %(x)].sketches['__profile__'].rectangle(point1=(length, height), 
        point2=(bottomcornerx, height-top_thickness))
    #Partition bottom section
    mdb.models['Model-%d' %(x)].sketches['__profile__'].rectangle(point1=(bottomcornerx,bottomcornery ), 
        point2=(length, bottom_thickness))
    #partition command
    mdb.models['Model-%d' %(x)].parts['Part-1'].PartitionFaceBySketch(faces=
        mdb.models['Model-%d' %(x)].parts['Part-1'].faces.findAt(((length/2,height/2,
        0),),), sketch=mdb.models['Model-%d' %(x)].sketches['__profile__'])
    del mdb.models['Model-%d' %(x)].sketches['__profile__']

    ##Materials creation (fixed)##
    #silicon
    mdb.models['Model-%d' %(x)].Material(name='Silicon')
    mdb.models['Model-%d' %(x)].materials['Silicon'].Elastic(table=((134E9, 
        0.22), ))
    #lithiated Area
    mdb.models['Model-%d' %(x)].Material(name='Lithiated Silicon')
    mdb.models['Model-%d' %(x)].materials['Lithiated Silicon'].Elastic(table=((
        41E9, 0.26), ))
    mdb.models['Model-%d' %(x)].materials['Lithiated Silicon'].Expansion(table=((
        2.61e-06, ), ))
    
    #Section creation  (fixed)##
    #section of lithium
    mdb.models['Model-%d' %(x)].HomogeneousSolidSection(material='Silicon', name=
        'Section-1', thickness=None)
    #section of lithiated area
    mdb.models['Model-%d' %(x)].HomogeneousSolidSection(material='Lithiated Silicon', 
        name='Section-2', thickness=None)
    #section set creation
    #top
    mdb.models['Model-%d' %(x)].parts['Part-1'].Set(faces=
        mdb.models['Model-%d' %(x)].parts['Part-1'].faces.findAt(((length/2,height-(top_thickness/2),
        0),)), name='topsectionset')
    #bottom
    mdb.models['Model-%d' %(x)].parts['Part-1'].Set(faces=
        mdb.models['Model-%d' %(x)].parts['Part-1'].faces.findAt(((length/2,bottom_thickness/2,
        0),)), name='bottomsectionset')
    #non-lithiated center
    mdb.models['Model-%d' %(x)].parts['Part-1'].Set(faces=
        mdb.models['Model-%d' %(x)].parts['Part-1'].faces.findAt(((length/2,height/2,
        0),)), name='midsectionset')
    #section assigment
    #assign silicon section 
    mdb.models['Model-%d' %(x)].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-%d' %(x)].parts['Part-1'].sets['midsectionset'], sectionName=
        'Section-1', thicknessAssignment=FROM_SECTION)

    #assign lithiated section
    mdb.models['Model-%d' %(x)].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-%d' %(x)].parts['Part-1'].sets['topsectionset'], sectionName=
        'Section-2', thicknessAssignment=FROM_SECTION)
    mdb.models['Model-%d' %(x)].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        mdb.models['Model-%d' %(x)].parts['Part-1'].sets['bottomsectionset'], sectionName=
        'Section-2', thicknessAssignment=FROM_SECTION)

    ##Create Nodal set along the reaction fronts##
    mdb.models['Model-%d' %(x)].parts['Part-1'].Set(edges=
        mdb.models['Model-%d' %(x)].parts['Part-1'].edges.getSequenceFromMask(('[#10 ]', 
        ), ), name='topreactionfront')
    mdb.models['Model-%d' %(x)].parts['Part-1'].Set(edges=
        mdb.models['Model-%d' %(x)].parts['Part-1'].edges.getSequenceFromMask(('[#1 ]', 
        ), ), name='bottomreactionfront')
    
    ##Assembly##
    #Create on instance
    mdb.models['Model-%d' %(x)].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-%d' %(x)].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-%d' %(x)].parts['Part-1'])
    #Add reference points for coupling
    distance_from_edges = 5
    #RP creation
    mdb.models['Model-%d' %(x)].rootAssembly.ReferencePoint(point=(
    0-distance_from_edges, height/2, 0.0))
    mdb.models['Model-%d' %(x)].rootAssembly.ReferencePoint(point=(
    length+ distance_from_edges, height/2, 0.0))
    #RP assignment
    RPleft_id = mdb.models['Model-%d' %(x)].rootAssembly.ReferencePoint(point=(0-distance_from_edges, height/2, 0.0)).id
    RPright_id = mdb.models['Model-%d' %(x)].rootAssembly.ReferencePoint(point=(length+ distance_from_edges, height/2, 0.0)).id

    
    
    ##Step Creation##
    #Buckling Step
    mdb.models['Model-%d' %(x)].BuckleStep(name='ThermalBuckling', numEigen=1, previous=
        'Initial', vectors=2)
    #Predfine Field
    temp = 19230
    mdb.models['Model-%d' %(x)].Temperature(createStepName='ThermalBuckling', 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
        UNIFORM, magnitudes=(temp, ), name='Predefined Field-top', region=Region(
        faces=mdb.models['Model-%d' %(x)].rootAssembly.instances['Part-1-1'].faces.findAt(((length/2,height-(top_thickness/2),
        0),))))
    mdb.models['Model-%d' %(x)].Temperature(createStepName='ThermalBuckling', 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, distributionType=
        UNIFORM, magnitudes=(temp, ), name='Predefined Field-bottom', region=Region(
        faces=mdb.models['Model-%d' %(x)].rootAssembly.instances['Part-1-1'].faces.findAt(((length/2,bottom_thickness/2,
        0),))))
    #Boundary Condition
    #In this case fixed both reference points 
    
    mdb.models['Model-%d' %(x)].DisplacementBC(amplitude=UNSET, 
    buckleCase=PERTURBATION_AND_BUCKLING, createStepName='ThermalBuckling', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'RP1', region=Region(referencePoints=(
    mdb.models['Model-%d' %(x)].rootAssembly.referencePoints[RPleft_id], ))
    , u1=0.0, u2=0.0, ur3=UNSET)
    
    mdb.models['Model-%d' %(x)].DisplacementBC(amplitude=UNSET, 
    buckleCase=PERTURBATION_AND_BUCKLING, createStepName='ThermalBuckling', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'RP2', region=Region(referencePoints=(
    mdb.models['Model-%d' %(x)].rootAssembly.referencePoints[RPright_id], ))
    , u1=0.0, u2=0.0, ur3=UNSET)

    #Coupling Constraints
    #left
    mdb.models['Model-%d' %(x)].Coupling(controlPoint=Region(
    referencePoints=(
    mdb.models['Model-%d' %(x)].rootAssembly.referencePoints[RPleft_id], ))
    , couplingType=DISTRIBUTING, influenceRadius=WHOLE_SURFACE, localCsys=None, 
    name='Constraint-1', surface=Region(
   edges=mdb.models['Model-%d' %(x)].rootAssembly.instances['Part-1-1'].edges.getByBoundingBox(0,-0.01,0,0,height,0)), u1=ON, u2=ON, ur3=ON, weightingMethod=UNIFORM)

    #right 
    mdb.models['Model-%d' %(x)].Coupling(controlPoint=Region(
    referencePoints=(
    mdb.models['Model-%d' %(x)].rootAssembly.referencePoints[RPright_id], ))
    , couplingType=DISTRIBUTING, influenceRadius=WHOLE_SURFACE, localCsys=None, 
    name='Constraint-2', surface=Region(
    edges=mdb.models['Model-%d' %(x)].rootAssembly.instances['Part-1-1'].edges.getByBoundingBox(length,-0.01,0,length,height,0)), u1=ON, u2=ON, ur3=ON, weightingMethod=UNIFORM)
    
    
    ##Creating a Mesh##
    #Element type
    mdb.models['Model-%d' %(x)].parts['Part-1'].setElementType(elemTypes=(ElemType(
        elemCode=CPE4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
        hourglassControl=DEFAULT, distortionControl=DEFAULT), ElemType(
        elemCode=CPE3, elemLibrary=STANDARD)), regions=(
        mdb.models['Model-%d' %(x)].parts['Part-1'].faces, ))
    #Seed Size (refer back to lithium size, as well as correlate it to the nodal set below)
    seedsize=30
    mdb.models['Model-%d' %(x)].parts['Part-1'].seedPart(deviationFactor=0.1, 
        minSizeFactor=0.1, size=seedsize)
    #Generate Mesh
    mdb.models['Model-%d' %(x)].parts['Part-1'].generateMesh()

    #Create Buckling Analysis Job
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-%d' %(x), modelPrint=OFF, name=
        'Job-%d' %(x), nodalOutputPrecision=SINGLE, queue=None, resultsFormat=ODB, 
        scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.models['Model-%d' %(x)].rootAssembly.regenerate()

    mdb.models['Model-%d' %(x)].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'E', 'U'))

    ##Create element set along the reaction fronts##
    #make sure 0.5 below is less than seedsize and correlated to deviation factor
    #bottom lithiated zone (gets bigger)
    mdb.models['Model-%d'%(x)].parts['Part-1'].Set(elements=
    mdb.models['Model-%d'%(x)].parts['Part-1'].elements.getByBoundingBox(-0.5,bottom_thickness-seedsize-5,0,length+5,bottom_thickness +1,0), name='bottomli')
    #bottom untransformed zone (only gets smaller)
    mdb.models['Model-%d'%(x)].parts['Part-1'].Set(elements=
    mdb.models['Model-%d'%(x)].parts['Part-1'].elements.getByBoundingBox(-0.5,bottom_thickness-0.01,0,length+5,bottom_thickness+seedsize+5,0), name='bottomsi')
    #top lithiated zone (gets bigger)
    mdb.models['Model-%d'%(x)].parts['Part-1'].Set(elements=
    mdb.models['Model-%d'%(x)].parts['Part-1'].elements.getByBoundingBox(-0.5,height-top_thickness-1,0,length+5,height-top_thickness+seedsize+5,0), name='topli')
    #top untransformed zone(only gets smaller)
    mdb.models['Model-%d'%(x)].parts['Part-1'].Set(elements=
    mdb.models['Model-%d'%(x)].parts['Part-1'].elements.getByBoundingBox(-0.5,height-top_thickness-seedsize-5,0,length+1,height-top_thickness+1,0), name='topsi')
    
    ##Submit Job-1##
    mdb.jobs['Job-%d' %(x)].submit()
    mdb.jobs['Job-%d' %(x)].waitForCompletion()

    ###Results###
    ##Data retrieval procedure

    ##Bottom Silicon##
    
    #Stress S
    #open ODB File
    odb = openOdb(path='Job-%d.odb' %(x))
    #Create a variable that refers to the - Step of the job
    selectstep = odb.steps['ThermalBuckling']
    # Create a variable that refers to the - Frame Number of the Step
    selectframe = selectstep.frames[1]
    # Create a variable that refers to the - the variable from the above frame of step - Stress 'S'
    selectvariable = selectframe.fieldOutputs['S']
    # Create a variable that refers to the element set 'bottomsi' (The set is associated with part instance 'PART-1-1'. USE CAPITAL for elementset name
    selectelement = odb.rootAssembly.instances['PART-1-1'].elementSets['BOTTOMSI']
    #Create a variable that refers to - reads the variable data corresponding to the element set region assigned (in this case stress for 'bottomsi')
    selectregion = selectvariable.getSubset(region=selectelement, position= INTEGRATION_POINT)
    #1. Print the field ouput data for all elements from the selected set

    s11bottomsi=[]
    s22bottomsi=[]
    s12bottomsi=[]
    txtlist=[]
    #This is just saving data per element to txt file for monitoring
    for s in selectregion.values:
        elementnumber = [s.elementLabel]
        elementS11 = [s.data[0]]
        elementS22 = [s.data[1]]
        elementS12 = [s.data[3]]
        elementposition =[s.position]
        #adding values to averaging list
        s11bottomsi.append(elementS11)
        s22bottomsi.append(elementS22)
        s12bottomsi.append(elementS12)
        #Printing them out to txt file (Trial using just braces/sets)      
        combineelementdata= [ 'Element Label = '+str(elementnumber), 'S11 = '+str(elementS11), 'S22 = '+str(elementS22), 'S12 = '+str(elementS12)]
        txtlist.extend(combineelementdata)
        np.savetxt('S'+'_elements_BottomSi_job%d.txt'%(x),txtlist,fmt='%s')
        
    #Calculate average of bottomsi of model
    sbottomsifinalaverage11 = sum(s11bottomsi)/len(s11bottomsi)
    sbottomsifinalaverage22 = sum(s22bottomsi)/len(s22bottomsi)
    sbottomsifinalaverage12 = sum(s12bottomsi)/len(s12bottomsi)
    bottomvalue = [sbottomsifinalaverage11, sbottomsifinalaverage22, sbottomsifinalaverage12]

    ##write averages to txt, weird occurances where list bottomvalue[m] not needed, only m when above i needed[] for each indication## 
    with open('C:/temp/S_bottom_SI-job-%d.txt' %(x),'w') as file:
            for m in bottomvalue:
                file.write('%s \n' %m)
                file.write('\n')
    #Strain E
    # Create a variable that refers to the - the variable from the above frame of step - Strain 'E'
    selectvariable = selectframe.fieldOutputs['E']
    # Create a variable that refers to the element set 'bottomsi' (The set is associated with part instance 'PART-1-1'. USE CAPITAL for elementset name
    selectelement = odb.rootAssembly.instances['PART-1-1'].elementSets['BOTTOMSI']
    #Create a variable that refers to - reads the variable data corresponding to the element set region assigned (in this case strain for 'bottomsi')
    selectregion = selectvariable.getSubset(region=selectelement, position= INTEGRATION_POINT)
    #1. Print the field ouput data for all elements from the selected set

    e11bottomsi=[]
    e22bottomsi=[]
    e12bottomsi=[]
    txtlist=[]
    #This is just saving data per element to txt file for monitoring
    for s in selectregion.values:
        elementnumber = [s.elementLabel]
        elementE11 = [s.data[0]]
        elementE22 = [s.data[1]]
        elementE12 = [s.data[3]]
        elementposition =[s.position]
        #adding values to averaging list
        e11bottomsi.append(elementE11)
        e22bottomsi.append(elementE22)
        e12bottomsi.append(elementE12)
        #Printing them out to txt file (Trial using just braces/sets)      
        combineelementdata= [ 'Element Label = '+str(elementnumber), 'E11 = '+str(elementE11), 'E22 = '+str(elementE22), 'E12 = '+str(elementE12)]
        txtlist.extend(combineelementdata)
        np.savetxt('E'+'_elements_BottomSi_job-%d.txt'%(x),txtlist,fmt='%s')
        
    #Calculate average of bottomsi of model
    Ebottomsifinalaverage11 = sum(e11bottomsi)/len(e11bottomsi)
    Ebottomsifinalaverage22 = sum(e22bottomsi)/len(e22bottomsi)
    Ebottomsifinalaverage12 = sum(e12bottomsi)/len(e12bottomsi)
    bottomvalue = [Ebottomsifinalaverage11, Ebottomsifinalaverage22, Ebottomsifinalaverage12]

    ##write averages to txt, weird occurances where list bottomvalue[m] not needed, only m when above i needed[] for each indication## 
    with open('C:/temp/E_bottom_SI-job%d.txt' %(x),'w') as file:
            for m in bottomvalue:
                file.write('%s \n' %m)
                file.write('\n')
                
                
    ##Bottom Lithiated Material##
                
    #Stress S
    # Create a variable that refers to the - the variable from the above frame of step - Stress 'S'
    selectvariable = selectframe.fieldOutputs['S']
    # Create a variable that refers to the element set 'bottomli' (The set is associated with part instance 'PART-1-1'. USE CAPITAL for elementset name
    selectelement = odb.rootAssembly.instances['PART-1-1'].elementSets['BOTTOMLI']
    #Create a variable that refers to - reads the variable data corresponding to the element set region assigned (in this case stress for 'bottomli')
    selectregion = selectvariable.getSubset(region=selectelement, position= INTEGRATION_POINT)
    #1. Print the field ouput data for all elements from the selected set

    s11bottomli=[]
    s22bottomli=[]
    s12bottomli=[]
    txtlist=[]
    #This is just saving data per element to txt file for monitoring
    for s in selectregion.values:
        elementnumber = [s.elementLabel]
        elementS11 = [s.data[0]]
        elementS22 = [s.data[1]]
        elementS12 = [s.data[3]]
        elementposition =[s.position]
        #adding values to averaging list
        s11bottomli.append(elementS11)
        s22bottomli.append(elementS22)
        s12bottomli.append(elementS12)
        #Printing them out to txt file (Trial using just braces/sets)      
        combineelementdata= [ 'Element Label = '+str(elementnumber), 'S11 = '+str(elementS11), 'S22 = '+str(elementS22), 'S12 = '+str(elementS12)]
        txtlist.extend(combineelementdata)
        np.savetxt('S'+'_elements_BottomLi_job%d.txt'%(x),txtlist,fmt='%s')
        
    #Calculate average of bottomli of model
    sbottomlifinalaverage11 = sum(s11bottomli)/len(s11bottomli)
    sbottomlifinalaverage22 = sum(s22bottomli)/len(s22bottomli)
    sbottomlifinalaverage12 = sum(s12bottomli)/len(s12bottomli)
    bottomvalue = [sbottomlifinalaverage11, sbottomlifinalaverage22, sbottomlifinalaverage12]

    ##write averages to txt, weird occurances where list bottomvalue[m] not needed, only m when above i needed[] for each indication## 
    with open('C:/temp/S_bottom_LI-job-%d.txt' %(x),'w') as file:
            for m in bottomvalue:
                file.write('%s \n' %m)
                file.write('\n')
                
    #Strain E
    # Create a variable that refers to the - the variable from the above frame of step - Strain 'E'
    selectvariable = selectframe.fieldOutputs['E']
    # Create a variable that refers to the element set 'bottomli' (The set is associated with part instance 'PART-1-1'. USE CAPITAL for elementset name
    selectelement = odb.rootAssembly.instances['PART-1-1'].elementSets['BOTTOMLI']
    #Create a variable that refers to - reads the variable data corresponding to the element set region assigned (in this case strain for 'bottomli')
    selectregion = selectvariable.getSubset(region=selectelement, position= INTEGRATION_POINT)
    #1. Print the field ouput data for all elements from the selected set

    e11bottomli=[]
    e22bottomli=[]
    e12bottomli=[]
    txtlist=[]
    #This is just saving data per element to txt file for monitoring
    for s in selectregion.values:
        elementnumber = [s.elementLabel]
        elementE11 = [s.data[0]]
        elementE22 = [s.data[1]]
        elementE12 = [s.data[3]]
        elementposition =[s.position]
        #adding values to averaging list
        e11bottomli.append(elementE11)
        e22bottomli.append(elementE22)
        e12bottomli.append(elementE12)
        #Printing them out to txt file (Trial using just braces/sets)      
        combineelementdata= [ 'Element Label = '+str(elementnumber), 'E11 = '+str(elementE11), 'E22 = '+str(elementE22), 'E12 = '+str(elementE12)]
        txtlist.extend(combineelementdata)
        np.savetxt('E'+'_elements_BottomLi_job-%d.txt'%(x),txtlist,fmt='%s')
        
    #Calculate average of bottomli of model
    Ebottomlifinalaverage11 = sum(e11bottomli)/len(e11bottomli)
    Ebottomlifinalaverage22 = sum(e22bottomli)/len(e22bottomli)
    Ebottomlifinalaverage12 = sum(e12bottomli)/len(e12bottomli)
    bottomvalue = [Ebottomlifinalaverage11, Ebottomlifinalaverage22, Ebottomlifinalaverage12]

    ##write averages to txt, weird occurances where list bottomvalue[m] not needed, only m when above i needed[] for each indication## 
    with open('C:/temp/E_bottom_LI-job%d.txt' %(x),'w') as file:
            for m in bottomvalue:
                file.write('%s \n' %m)
                file.write('\n')
    

    ##Top Silicon##
                
    #Stress S
    # Create a variable that refers to the - the variable from the above frame of step - Stress 'S'
    selectvariable = selectframe.fieldOutputs['S']
    # Create a variable that refers to the element set 'TopSi' (The set is associated with part instance 'PART-1-1'. USE CAPITAL for elementset name
    selectelement = odb.rootAssembly.instances['PART-1-1'].elementSets['TOPSI']
    #Create a variable that refers to - reads the variable data corresponding to the element set region assigned (in this case stress for 'topsi')
    selectregion = selectvariable.getSubset(region=selectelement, position= INTEGRATION_POINT)
    #1. Print the field ouput data for all elements from the selected set

    s11topsi=[]
    s22topsi=[]
    s12topsi=[]
    txtlist=[]
    #This is just saving data per element to txt file for monitoring
    for s in selectregion.values:
        elementnumber = [s.elementLabel]
        elementS11 = [s.data[0]]
        elementS22 = [s.data[1]]
        elementS12 = [s.data[3]]
        elementposition =[s.position]
        #adding values to averaging list
        s11topsi.append(elementS11)
        s22topsi.append(elementS22)
        s12topsi.append(elementS12)
        #Printing them out to txt file (Trial using just braces/sets)      
        combineelementdata= [ 'Element Label = '+str(elementnumber), 'S11 = '+str(elementS11), 'S22 = '+str(elementS22), 'S12 = '+str(elementS12)]
        txtlist.extend(combineelementdata)
        np.savetxt('S'+'_elements_TOPSI_job%d.txt'%(x),txtlist,fmt='%s')
        
    #Calculate average of topsi of model
    stopsifinalaverage11 = sum(s11topsi)/len(s11topsi)
    stopsifinalaverage22 = sum(s22topsi)/len(s22topsi)
    stopsifinalaverage12 = sum(s12topsi)/len(s12topsi)
    bottomvalue = [stopsifinalaverage11, stopsifinalaverage22, stopsifinalaverage12]

    ##write averages to txt, weird occurances where list bottomvalue[m] not needed, only m when above i needed[] for each indication## 
    with open('C:/temp/S_top_SI-job-%d.txt' %(x),'w') as file:
            for m in bottomvalue:
                file.write('%s \n' %m)
                file.write('\n')
                
    #Strain E
    # Create a variable that refers to the - the variable from the above frame of step - Strain 'E'
    selectvariable = selectframe.fieldOutputs['E']
    # Create a variable that refers to the element set 'topsi' (The set is associated with part instance 'PART-1-1'. USE CAPITAL for elementset name
    selectelement = odb.rootAssembly.instances['PART-1-1'].elementSets['TOPSI']
    #Create a variable that refers to - reads the variable data corresponding to the element set region assigned (in this case strain for 'topsi')
    selectregion = selectvariable.getSubset(region=selectelement, position= INTEGRATION_POINT)
    #1. Print the field ouput data for all elements from the selected set

    e11topsi=[]
    e22topsi=[]
    e12topsi=[]
    txtlist=[]
    #This is just saving data per element to txt file for monitoring
    for s in selectregion.values:
        elementnumber = [s.elementLabel]
        elementE11 = [s.data[0]]
        elementE22 = [s.data[1]]
        elementE12 = [s.data[3]]
        elementposition =[s.position]
        #adding values to averaging list
        e11topsi.append(elementE11)
        e22topsi.append(elementE22)
        e12topsi.append(elementE12)
        #Printing them out to txt file (Trial using just braces/sets)      
        combineelementdata= [ 'Element Label = '+str(elementnumber), 'E11 = '+str(elementE11), 'E22 = '+str(elementE22), 'E12 = '+str(elementE12)]
        txtlist.extend(combineelementdata)
        np.savetxt('E'+'_elements_TopSi_job-%d.txt'%(x),txtlist,fmt='%s')
        
    #Calculate average of bottomsi of model
    Etopsifinalaverage11 = sum(e11topsi)/len(e11topsi)
    Etopsifinalaverage22 = sum(e22topsi)/len(e22topsi)
    Etopsifinalaverage12 = sum(e12topsi)/len(e12topsi)
    bottomvalue = [Etopsifinalaverage11, Etopsifinalaverage22, Etopsifinalaverage12]

    ##write averages to txt, weird occurances where list bottomvalue[m] not needed, only m when above i needed[] for each indication## 
    with open('C:/temp/E_top_SI-job%d.txt' %(x),'w') as file:
            for m in bottomvalue:
                file.write('%s \n' %m)
                file.write('\n')

    ##Top Lithiated Material##
                
    #Stress S
    # Create a variable that refers to the - the variable from the above frame of step - Stress 'S'
    selectvariable = selectframe.fieldOutputs['S']
    # Create a variable that refers to the element set 'TopLi' (The set is associated with part instance 'PART-1-1'. USE CAPITAL for elementset name
    selectelement = odb.rootAssembly.instances['PART-1-1'].elementSets['TOPLI']
    #Create a variable that refers to - reads the variable data corresponding to the element set region assigned (in this case stress for 'topli')
    selectregion = selectvariable.getSubset(region=selectelement, position= INTEGRATION_POINT)
    #1. Print the field ouput data for all elements from the selected set

    s11topli=[]
    s22topli=[]
    s12topli=[]
    txtlist=[]
    #This is just saving data per element to txt file for monitoring
    for s in selectregion.values:
        elementnumber = [s.elementLabel]
        elementS11 = [s.data[0]]
        elementS22 = [s.data[1]]
        elementS12 = [s.data[3]]
        elementposition =[s.position]
        #adding values to averaging list
        s11topli.append(elementS11)
        s22topli.append(elementS22)
        s12topli.append(elementS12)
        #Printing them out to txt file (Trial using just braces/sets)      
        combineelementdata= [ 'Element Label = '+str(elementnumber), 'S11 = '+str(elementS11), 'S22 = '+str(elementS22), 'S12 = '+str(elementS12)]
        txtlist.extend(combineelementdata)
        np.savetxt('S'+'_elements_TOPLI_job%d.txt'%(x),txtlist,fmt='%s')
        
    #Calculate average of topsi of model
    stoplifinalaverage11 = sum(s11topli)/len(s11topli)
    stoplifinalaverage22 = sum(s22topli)/len(s22topli)
    stoplifinalaverage12 = sum(s12topli)/len(s12topli)
    bottomvalue = [stoplifinalaverage11, stoplifinalaverage22, stoplifinalaverage12]

    ##write averages to txt, weird occurances where list bottomvalue[m] not needed, only m when above i needed[] for each indication## 
    with open('C:/temp/S_top_LI-job-%d.txt' %(x),'w') as file:
            for m in bottomvalue:
                file.write('%s \n' %m)
                file.write('\n')
                
    #Strain E
    # Create a variable that refers to the - the variable from the above frame of step - Strain 'E'
    selectvariable = selectframe.fieldOutputs['E']
    # Create a variable that refers to the element set 'topli' (The set is associated with part instance 'PART-1-1'. USE CAPITAL for elementset name
    selectelement = odb.rootAssembly.instances['PART-1-1'].elementSets['TOPLI']
    #Create a variable that refers to - reads the variable data corresponding to the element set region assigned (in this case strain for 'topli')
    selectregion = selectvariable.getSubset(region=selectelement, position= INTEGRATION_POINT)
    #1. Print the field ouput data for all elements from the selected set

    e11topli=[]
    e22topli=[]
    e12topli=[]
    txtlist=[]
    #This is just saving data per element to txt file for monitoring
    for s in selectregion.values:
        elementnumber = [s.elementLabel]
        elementE11 = [s.data[0]]
        elementE22 = [s.data[1]]
        elementE12 = [s.data[3]]
        elementposition =[s.position]
        #adding values to averaging list
        e11topli.append(elementE11)
        e22topli.append(elementE22)
        e12topli.append(elementE12)
        #Printing them out to txt file (Trial using just braces/sets)      
        combineelementdata= [ 'Element Label = '+str(elementnumber), 'E11 = '+str(elementE11), 'E22 = '+str(elementE22), 'E12 = '+str(elementE12)]
        txtlist.extend(combineelementdata)
        np.savetxt('E'+'_elements_TopLi_job-%d.txt'%(x),txtlist,fmt='%s')
        
    #Calculate average of bottomsi of model
    Etoplifinalaverage11 = sum(e11topli)/len(e11topli)
    Etoplifinalaverage22 = sum(e22topli)/len(e22topli)
    Etoplifinalaverage12 = sum(e12topli)/len(e12topli)
    bottomvalue = [Etoplifinalaverage11, Etoplifinalaverage22, Etoplifinalaverage12]

    ##write averages to txt, weird occurances where list bottomvalue[m] not needed, only m when above i needed[] for each indication## 
    with open('C:/temp/E_top_LI-job%d.txt' %(x),'w') as file:
            for m in bottomvalue:
                file.write('%s \n' %m)
    
    
    ###calculate new reaction height (testing with same thickness on both sides)###
    #Variables
    #(1)kinetic coefficient(reaction front coefficient) m/s
    kstar= 8.6e-8
    #(2)reference concentration of the diffusing constituent taken as a solubility mol/m^3
    cstar= 53000
    #(3)stoichiometric coefficients of initial solid material
    nminus = 4.0/15
    #(4)Molar Mass of initial solid material kg/Mol
    mminus = 0.0280855
    #(5)reference mass density kg/m3
    rhominus = 2285
    #combined variables (3),(4),(5)
    initialvariable =(nminus*mminus)/rhominus
    #chemical transformation strain (Thermal in this case)
    ech= 0.05
    #energy parameter defined by the chemical energies of stress-free solid constituent J/m3
    bigy= 5e9

    #Bottom new reaction height#
    #bulk terms to fit into the equation 
    bottom_untransformedrho_untransformedstrain = sbottomsifinalaverage11*Ebottomsifinalaverage11+2*(sbottomsifinalaverage12*Ebottomsifinalaverage12)+sbottomsifinalaverage22*Ebottomsifinalaverage22

    bottom_transformedrho_transformedstrainwithchem = sbottomlifinalaverage11*(Ebottomlifinalaverage11-ech)+2*sbottomlifinalaverage12*Ebottomlifinalaverage12+sbottomlifinalaverage22*(Ebottomlifinalaverage22-ech)

    bottom_untransformedrho_transformedstrainuntransformedstrain = sbottomsifinalaverage11*(Ebottomlifinalaverage11-Ebottomsifinalaverage11)+2*sbottomsifinalaverage12*(Ebottomlifinalaverage12-Ebottomsifinalaverage12)+sbottomsifinalaverage22*(Ebottomlifinalaverage22-Ebottomsifinalaverage22)

    #Contributions of mechanical and chemical energies
    bottombigx = bigy +0.5*(bottom_untransformedrho_untransformedstrain)-0.5*(bottom_transformedrho_transformedstrainwithchem)+(bottom_untransformedrho_transformedstrainuntransformedstrain)

    #normal component of the chemical affinity tensor within a small strain approximation:
    ANN = initialvariable*bottombigx
    #Temperature in K (obtained from E = aT)
    temp= 293
    RT=8.314*temp
    expcomponent= -(ANN/RT)
    EXP= math.exp(expcomponent)
    #Reaction rate reformulated from kinetic equation by Glansdorff
    bottomwn = kstar*cstar*(1-EXP)
    #Reaction front velocity vwn cm/s
    bottom_vwn_m = initialvariable*bottomwn
    #Reaction front velocity in nm/s
    bottom_vwn = bottom_vwn_m/1E-9
    #distance travelled in this time step:
    bottomdistancetravel = bottom_vwn*timestep
    #new bottom recation front height m 
    new_bottom_reaction_height = bottomdistancetravel +bottom_thickness
    print("New bottom reaction front velocity is %.9f nm/s" %(bottom_vwn))
    print("New bottom reaction height after %.2f seconds: %.9f nm" %(timestep,new_bottom_reaction_height))

    #Top new reaction height#
    #bulk terms to fit into the equation 
    top_untransformedrho_untransformedstrain = stopsifinalaverage11*Etopsifinalaverage11+2*(stopsifinalaverage12*Etopsifinalaverage12)+stopsifinalaverage22*Etopsifinalaverage22

    top_transformedrho_transformedstrainwithchem = stoplifinalaverage11*(Etoplifinalaverage11-ech)+2*stoplifinalaverage12*Etoplifinalaverage12+stoplifinalaverage22*(Etoplifinalaverage22-ech)

    top_untransformedrho_transformedstrainuntransformedstrain = stopsifinalaverage11*(Etoplifinalaverage11-Etopsifinalaverage11)+2*stopsifinalaverage12*(Etoplifinalaverage12-Etopsifinalaverage12)+stopsifinalaverage22*(Etoplifinalaverage22-Etopsifinalaverage22)

    #Contributions of mechanical and chemical energies
    topbigx = bigy +0.5*(top_untransformedrho_untransformedstrain)-0.5*(top_transformedrho_transformedstrainwithchem)+(top_untransformedrho_transformedstrainuntransformedstrain)

    #normal component of the chemical affinity tensor within a small strain approximation:
    ANN = initialvariable*topbigx
    #Temperature in K (obtained from E = aT)
    expcomponent= -(ANN/RT)
    EXP= math.exp(expcomponent)
    #Reaction rate reformulated from kinetic equation by Glansdorff
    topwn = kstar*cstar*(1-EXP)
    #Reaction front velocity vwn m/s
    top_vwn_m = initialvariable*topwn
    #Reaction front velocity nm/s
    top_vwn=top_vwn_m/1E-9
    #distance travelled in this time step:
    topdistancetravel = top_vwn*timestep
    #new top recation front height m 
    new_top_reaction_height = height-top_thickness-topdistancetravel
    print("New top reaction front velocity is %.9f nm/s" %(top_vwn))
    print"New top reaction height after %.2f seconds: %.9f nm" %(timestep,new_top_reaction_height)
    

    #eigenvalue check if it buckled
    #open ODB File
    odb = openOdb(path='Job-%d.odb' %(x))
    #Create a variable that refers to the - Step of the job
    selectstep = odb.steps['ThermalBuckling']
    # Create a variable that refers to the - Frame Number of the Step
    selectframe = selectstep.frames[1]
    # Create a variable that refers to the - the variable from the above frame of step - eigenvalue
    eigenvalueandmode = selectframe.description
    eigenvaluestr = eigenvalueandmode[28:48]
    eigenvalue=float(eigenvaluestr)
    print("Eigenvalue of this time step is :%.9f" %eigenvalue)
    
    #Relative thicknesses
    toprelativethickness = top_thickness/(height/2)
    bottomrelativethickness =  bottom_thickness/(height/2)
    
    #new thicknesses
    top_thickness = top_thickness+topdistancetravel
    bottom_thickness = bottom_thickness+ bottomdistancetravel

    #Exporting the data for plotting
    topreactionvtxtlist.append(top_vwn)
    bottomreactionvtxtlist.append(bottom_vwn)
    topreactionrelativethicknesstxtlist.append(toprelativethickness)
    bottomreactionrelativethicknesstxtlist.append(bottomrelativethickness)
    

else:
    print('critical thickness of transformed layer reached')
    list_tuples = list(zip(topreactionvtxtlist, topreactionrelativethicknesstxtlist, bottomreactionvtxtlist, bottomreactionrelativethicknesstxtlist))
    arr = np.asarray (list_tuples)
    print ("Array:", arr)
    np.savetxt("testing.csv", arr, delimiter=",")





#check week 20 for notes 
    
