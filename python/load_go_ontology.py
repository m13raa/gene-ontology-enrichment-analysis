#Create a GO DAG object from the OBO file, containing the GO ontology
from goatools.obo_parser import GODag
obo_file = "data/raw/go-basic.obo"
go_dag = GODag(obo_file)
#Checks
print("GO ontology loaded")
print(len(go_dag))
#Load ontology
go_dag = GODag(str(obo_file))
print(F"\nLoaded {len(go_dag):,} GO terms.\n")
#Example of a GO term
go_term = go_dag["GO:0006915"]
print("GO ID:", go_term.id)
print("Name:", go_term.name)
print("Namespace:", go_term.namespace)
print("\nOntology details/info")
print("\nParents:")
print("Number of parents:", len(go_term.parents))
for parent in go_term.parents:
    print(f"  - {parent.id}: {parent.name}")
print("\nChildren:")
print("Number of children:", len(go_term.children))
for child in go_term.children:
    print(f"  - {child.id}: {child.name}")
print("Depth:", go_term.depth)
print("Level:", go_term.level)