from paraview import simple
reader = simple.OpenDataFile("/phreeqc/porous_compositional_1.vtu")
writer = simple.CreateWriter("/phreeqc/porous_compositional_1.csv", reader)
writer.WriteAllTimeSteps = 1
writer.FieldAssociation = "Points"
writer.UpdatePipeline()
del writer
