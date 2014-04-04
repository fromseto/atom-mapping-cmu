#include <iostream>
#include <vector>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/reaction.h>

using namespace std;
using namespace OpenBabel;

/* global definition */
/* reactants*/
std::vector<int> rtypes;
std::vector<int> rbonds;
std::vector<int> rtetra;
std::vector<int> rta;
std::vector<int> rtt;
std::vector<int> rdb;
std::vector<int> rda;

/* products*/
std::vector<int> ptypes;
std::vector<int> pbonds;
std::vector<int> ptetra;
std::vector<int> pta;
std::vector<int> ptt;
std::vector<int> pdb;
std::vector<int> pda;

// subfunction prototype
int get_info_of_each_reactant(shared_ptr<OBMol> m);
int get_info_of_each_product(shared_ptr<OBMol> m);
void get_neighbor_reactant(int num_atoms);
void get_neighbor_product(int num_atoms);
//void print_to_file(ofstream myfile);



int main(int argc,char **argv)
{

  //OpenBabel::OBConversion conv;
  OBConversion conv;
  //conv.SetInFormat("sdf");
  conv.SetInFormat("rxn");
  //OpenBabel::OBMol mol;
  //OBMol mol;
  OBReaction rxn;
  
  
  /*
  bool notatend = conv.ReadFile(&mol, "my_first.sdf");
  
  while (notatend){
    std::cout << "Molecular weight: "<< mol.GetMolWt() << std::endl;
    mol.Clear();
    notatend = conv.Read(&mol);
    }
    */
    //bool notatend = conv.ReadFile(&rxn, "my_first.rxn");
    bool notatend = conv.ReadFile(&rxn, argv[1]);
    ofstream myfile ("gams_input_example.txt");

    while (notatend){
        int num_reactants = rxn.NumReactants();
        int num_products = rxn.NumProducts();
        std::cout << "Number of reactants: "<< num_reactants << std::endl;
        std::cout << "Number of Productss: "<< num_products << std::endl;
        shared_ptr<OBMol> m = rxn.GetReactant(0);
        shared_ptr<OBMol> m1 = rxn.GetReactant(1);

        shared_ptr<OBMol> n = rxn.GetProduct(0);
        shared_ptr<OBMol> n1 = rxn.GetProduct(1);

        std::cout << "Number of atoms "<< (*m).NumAtoms() << std::endl;
        std::cout << "Formula "<< (*m).GetFormula() << std::endl;
        std::cout << "Number of atoms "<< (*m1).NumAtoms() << std::endl;
        std::cout << "Formula "<< (*m1).GetFormula() << std::endl;

        std::cout << "Number of atoms "<< (*n).NumAtoms() << std::endl;
        std::cout << "Formula "<< (*n).GetFormula() << std::endl;

        std::cout << "Number of atoms "<< (*n1).NumAtoms() << std::endl;
        std::cout << "Formula "<< (*n1).GetFormula() << std::endl;

        int num1, num2;
        num1 = get_info_of_each_reactant(m);
        num1 = get_info_of_each_reactant(m1);
        get_neighbor_reactant(num1);

        num2 = get_info_of_each_product(n);
        num2 = get_info_of_each_product(n1);
        get_neighbor_product(num2);

        //print_to_file(myfile);
        
    // out put number of atoms
        myfile << "Set i /a1 * a"<<(*m).NumAtoms() + (*m1).NumAtoms()<<"/;\n";
        myfile << "Alias(i,j,k,l,m,n,o,p);\n";

        //myfile << "* reactants\n";

    // output all atom types
        myfile << "Parameter rtypes(i) 'types of atoms'/";
        if (rtypes.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rtypes.size(); ix++)
            {
                if (ix==rtypes.size() - 1)
                {
                    myfile << "a"<< ix+1 <<" "<<rtypes[ix]<<"/;\n";
                } else {
                    myfile << "a"<< ix+1 <<" "<<rtypes[ix]<<", ";
                }   
            }
        } else {
            myfile << " /;\n";
        }

        myfile << "Parameter ptypes(k) 'types of atoms'/";
        if (ptypes.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != ptypes.size(); ix++)
            {
                if (ix==ptypes.size() - 1)
                {
                    myfile << "a"<< ix+1 <<" "<<ptypes[ix]<<"/;\n\n";
                } else {
                    myfile << "a"<< ix+1 <<" "<<ptypes[ix]<<", ";
                }   
            }
        } else {
            myfile << " /;\n\n";
        }

    // output all bond connection between two atoms for gams file
        myfile << "Set rbonds(i,j) 'bonds btween two atoms ' /";

        if (rbonds.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rbonds.size(); ix += 2)
            {
                if (ix==rbonds.size() - 2)
                {
                    myfile << "a"<<rbonds[ix]<<"."<< "a"<<rbonds[ix+1]<<"/;\n";
                } else {
                    myfile << "a"<<rbonds[ix]<<"."<< "a"<<rbonds[ix+1]<<", ";
                }   
            }
        } else {
            myfile << " /;\n";
        }

        myfile << "Set pbonds(k,l) 'bonds btween two atoms' /";

        if (pbonds.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != pbonds.size(); ix += 2)
            {
                if (ix==rbonds.size() - 2)
                {
                    myfile << "a"<<pbonds[ix]<<"."<< "a"<<pbonds[ix+1]<<"/;\n\n";
                } else {
                    myfile << "a"<<pbonds[ix]<<"."<< "a"<<pbonds[ix+1]<<", ";
                }   
            }
        } else {
            myfile << " /;\n\n";
        }

// output tetra atoms
        myfile << "Set adjr /t1 * t4/;\n";
        myfile << "Alias(adjr,adjp);\n\n";
        myfile << "Set rtetra(i) 'index of tetra atom' /";
        if (rtetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rtetra.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< rtetra[ix] <<" /;\n";
                } else {
                    myfile << "a"<< rtetra[ix] <<" , ";
                }   
            }
        } else {
            myfile << " /;\n";
        }

        myfile << "Set ptetra(k) 'index of tetra atom' /";
        if (rtetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != ptetra.size(); ix++)
            {
                if (ix==ptetra.size() - 1)
                {
                    myfile << "a"<< ptetra[ix] <<" /;\n\n";
                } else {
                    myfile << "a"<< ptetra[ix] <<" , ";
                }   
            }
        } else {
            myfile << " /;\n\n";
        }


// output four neighoring atoms of each tetra atoms
        myfile << "Parameter rta(i,adjr) 'neighoring atoms of each tetra atoms'/";
        if (rtetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rtetra.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< rtetra[ix] <<".t1 "<< rta[4*ix]<<" , "
                    << "a"<< rtetra[ix] <<".t2 "<< rta[4*ix+1]<<" , "
                    << "a"<< rtetra[ix] <<".t3 "<< rta[4*ix+2]<<" , "
                    << "a"<< rtetra[ix] <<".t4 "<< rta[4*ix+3]<<" /;\n";
                } else {
                    myfile << "a"<< rtetra[ix] <<".t1 "<< rta[4*ix]<<" , "
                    << "a"<< rtetra[ix] <<".t2 "<< rta[4*ix+1]<<" , "
                    << "a"<< rtetra[ix] <<".t3 "<< rta[4*ix+2]<<" , "
                    << "a"<< rtetra[ix] <<".t4 "<< rta[4*ix+3]<<" , ";
                }   
            }
        } else {
            myfile << " /;\n";
        }

        myfile << "Parameter pta(k,adjp) 'neighoring atoms of each tetra atoms'/";
        if (ptetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != ptetra.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< ptetra[ix] <<".t1 "<< pta[4*ix]<<" , "
                    << "a"<< ptetra[ix] <<".t2 "<< pta[4*ix+1]<<" , "
                    << "a"<< ptetra[ix] <<".t3 "<< pta[4*ix+2]<<" , "
                    << "a"<< ptetra[ix] <<".t4 "<< pta[4*ix+3]<<" /;\n\n";
                } else {
                    myfile << "a"<< ptetra[ix] <<".t1 "<< pta[4*ix]<<" , "
                    << "a"<< ptetra[ix] <<".t2 "<< pta[4*ix+1]<<" , "
                    << "a"<< ptetra[ix] <<".t3 "<< pta[4*ix+2]<<" , "
                    << "a"<< ptetra[ix] <<".t4 "<< pta[4*ix+3]<<" , ";
                }   
            }
        } else {
            myfile << " /;\n\n";
        }

// output types of tetra atoms for gams file
        myfile << "Parameter rtt(i) 'types of tetra atoms' /";
        if (rtetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rtetra.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< rtetra[ix] <<" 6/;\n";
                } else {
                    myfile << "a"<< rtetra[ix] <<" 6, ";
                }   
            }
        } else {
            myfile << " /;\n";
        }

        myfile << "Parameter ptt(k) 'types of tetra atoms' /";
        if (ptetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != ptetra.size(); ix++)
            {
                if (ix==ptetra.size() - 1)
                {
                    myfile << "a"<< ptetra[ix] <<" 6/;\n\n";
                } else {
                    myfile << "a"<< ptetra[ix] <<" 6, ";
                }   
            }
        } else {
            myfile << " /;\n\n";
        }

// output double bonds for gams file
        myfile << "Set rdb(i,j) 'double bonds between two atom index'/";
        if (rdb.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rdb.size(); ix += 2)
            {
                if (ix==rdb.size() - 2)
                {
                    myfile << "a"<<rdb[ix]<<"."<< "a"<<rdb[ix+1]<<"/;\n";
                }
                myfile << "a"<<rdb[ix]<<"."<< "a"<<rdb[ix+1]<<", ";   
            }
        } else {
            myfile << " /;\n";
        }

        myfile << "Set pdb(k,l) 'double bonds between two atom index'/";
        if (rdb.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != pdb.size(); ix += 2)
            {
                if (ix==pdb.size() - 2)
                {
                    myfile << "a"<<pdb[ix]<<"."<< "a"<<pdb[ix+1]<<"/;\n\n";
                }
                myfile << "a"<<pdb[ix]<<"."<< "a"<<pdb[ix+1]<<", ";   
            }
        } else {
            myfile << " /;\n\n";
        }


// output four neighoring atoms of double bonds atoms
        myfile << "Parameter rda(i,j,adjr) 'neighoring atoms of double bonds'/";
        if (rdb.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rdb.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< rdb[ix] <<".t1 "<< rta[4*ix]<<" , "
                    << "a"<< rdb[ix] <<".t2 "<< rta[4*ix+1]<<" /;\n";
                } else {
                    myfile << "a"<< rdb[ix] <<".t1 "<< rta[4*ix]<<" , "
                    << "a"<< rdb[ix] <<".t2 "<< rta[4*ix+1]<<" , ";
                }   
            }
        } else {
            myfile << " /;\n";
        } 

        myfile << "Parameter pda(k,l,adjp) 'neighoring atoms of double bonds'/";
        if (rdb.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != pdb.size(); ix++)
            {
                if (ix==ptetra.size() - 1)
                {
                    myfile << "a"<< pdb[ix] <<".t1 "<< pta[4*ix]<<" , "
                    << "a"<< pdb[ix] <<".t2 "<< pta[4*ix+1]<<" /;\n\n";
                } else {
                    myfile << "a"<< pdb[ix] <<".t1 "<< pta[4*ix]<<" , "
                    << "a"<< pdb[ix] <<".t2 "<< pta[4*ix+1]<<" , ";
                }   
            }
        } else {
            myfile << " /;\n\n";
        } 

        rxn.Clear();
        notatend = conv.Read(&rxn);
    }

    myfile.close();
    return 0;
}

int get_info_of_each_reactant(shared_ptr<OBMol> m)
{        // get information of atoms
        // atom index start with 1 in OpenBabel
    int num_atoms = (*m).NumAtoms();
    //std::vector<int> rtypes;
    int magic_number = rtypes.size();
    for (int index = 1; index <= num_atoms; index++)
    {
        OBAtom atom = *((*m).GetAtom(index));
        std::cout << "Type of atom "<< *(atom.GetType()) << atom.GetAtomicMass() << std::endl;
        rtypes.push_back(int(atom.GetAtomicNum()));

    }        

        // get information of bonds
        // bond index start with 0 in OpenBabel
    int num_bonds = (*m).NumBonds();
    std::cout << "Number of bonds "<< num_bonds << std::endl;
    //std::vector<int> rbonds;
    //std::vector<int> rdb;

    for (int index = 0; index < num_bonds; index++)
    {
        OBBond bond = *((*m).GetBond(index));
        std::cout << "bond atom"<< (*(bond.GetBeginAtom())).GetIdx()+magic_number 
        << *((*(bond.GetBeginAtom())).GetType()) << " bond atom"
        << (*(bond.GetEndAtom())).GetIdx() +magic_number << *((*(bond.GetEndAtom())).GetType())
        << std::endl;

            // bond connection between two atoms
        rbonds.push_back((*(bond.GetBeginAtom())).GetIdx()+magic_number);
        rbonds.push_back((*(bond.GetEndAtom())).GetIdx()+magic_number);


            // double bonds
        if (bond.IsDouble()) 
        {
            rdb.push_back((*(bond.GetBeginAtom())).GetIdx()+magic_number);
            rdb.push_back((*(bond.GetEndAtom())).GetIdx()+magic_number);
        }
    }

    return magic_number;
}

int get_info_of_each_product(shared_ptr<OBMol> m)
{        // get information of atoms
        // atom index start with 1 in OpenBabel
    int num_atoms = (*m).NumAtoms();
    //std::vector<int> rtypes;
    int magic_number = ptypes.size();
    for (int index = 1; index <= num_atoms; index++)
    {
        OBAtom atom = *((*m).GetAtom(index));
        std::cout << "Type of atom "<< *(atom.GetType()) << atom.GetAtomicMass() << std::endl;
        ptypes.push_back(int(atom.GetAtomicNum()));

    }        

        // get information of bonds
        // bond index start with 0 in OpenBabel
    int num_bonds = (*m).NumBonds();
    std::cout << "Number of bonds "<< num_bonds << std::endl;
    //std::vector<int> rbonds;
    //std::vector<int> rdb;

    for (int index = 0; index < num_bonds; index++)
    {
        OBBond bond = *((*m).GetBond(index));
        std::cout << "bond atom"<< (*(bond.GetBeginAtom())).GetIdx()+magic_number 
        << *((*(bond.GetBeginAtom())).GetType()) << " bond atom"
        << (*(bond.GetEndAtom())).GetIdx() +magic_number << *((*(bond.GetEndAtom())).GetType())
        << std::endl;

            // bond connection between two atoms
        pbonds.push_back((*(bond.GetBeginAtom())).GetIdx()+magic_number);
        pbonds.push_back((*(bond.GetEndAtom())).GetIdx()+magic_number);


            // double bonds
        if (bond.IsDouble()) 
        {
            pdb.push_back((*(bond.GetBeginAtom())).GetIdx()+magic_number);
            pdb.push_back((*(bond.GetEndAtom())).GetIdx()+magic_number);
        }
    }

    return magic_number;
}

void get_neighbor_reactant(int num_atoms)
{
        // determine which is tetra atom
    //std::vector<int> rtetra;
        //std::vector<int> rtt;
    for (int i = 0; i <= num_atoms ; ++i)
    {
        int tetra_count = 0;
        if (rtypes[i] == 6) 
        {
            for(vector<int>::size_type ix = 0; ix != rbonds.size(); ix ++){
                if (rbonds[ix] == (i + 1)) 
                {
                    tetra_count += 1;
                }
            }
            if (tetra_count == 4) {
               rtetra.push_back(i+1);
           } 
       }
   }

     // find four neighoring atoms to double bonds
   //std::vector<int> rda;
   for (vector<int>::size_type ix = 0; ix != rdb.size(); ix ++)
   {
    for(vector<int>::size_type iy = 0; iy != rbonds.size(); iy += 2){
        if (rdb[ix] == rbonds[iy] && rdb[ix+1] != rbonds[iy])
        {
            rda.push_back(rbonds[iy+1]);
        }
    }
    for(vector<int>::size_type iy = 1; iy != (rbonds.size() + 1); iy += 2){
        if (rdb[ix] == rbonds[iy] && rdb[ix+1] != rbonds[iy]) {
            rda.push_back(rbonds[iy-1]);
        }
    }
}

    // find four neighoring atoms to tetra atom
//std::vector<int> rta;
for (vector<int>::size_type ix = 0; ix != rtetra.size(); ix ++)
{
 for(vector<int>::size_type iy = 0; iy != rbonds.size(); iy += 2){
    if (rtetra[ix] == rbonds[iy])
    {
        rta.push_back(rbonds[iy+1]);
    }
}
for(vector<int>::size_type iy = 1; iy != (rbonds.size() + 1); iy += 2){
    if (rtetra[ix] == rbonds[iy])
    {
        rta.push_back(rbonds[iy-1]);
    }
}
}

}

void get_neighbor_product(int num_atoms)
{
        // determine which is tetra atom
    //std::vector<int> rtetra;
        //std::vector<int> rtt;
    for (int i = 0; i <= num_atoms ; ++i)
    {
        int tetra_count = 0;
        if (ptypes[i] == 6) 
        {
            for(vector<int>::size_type ix = 0; ix != pbonds.size(); ix ++){
                if (pbonds[ix] == (i + 1)) 
                {
                    tetra_count += 1;
                }
            }
            if (tetra_count == 4) {
               ptetra.push_back(i+1);
           } 
       }
   }

     // find four neighoring atoms to double bonds
   //std::vector<int> rda;
   for (vector<int>::size_type ix = 0; ix != pdb.size(); ix ++)
   {
    for(vector<int>::size_type iy = 0; iy != pbonds.size(); iy += 2){
        if (pdb[ix] == pbonds[iy] && pdb[ix+1] != pbonds[iy])
        {
            pda.push_back(pbonds[iy+1]);
        }
    }
    for(vector<int>::size_type iy = 1; iy != (pbonds.size() + 1); iy += 2){
        if (pdb[ix] == pbonds[iy] && pdb[ix+1] != pbonds[iy]) {
            pda.push_back(pbonds[iy-1]);
        }
    }
}

    // find four neighoring atoms to tetra atom
//std::vector<int> rta;
for (vector<int>::size_type ix = 0; ix != ptetra.size(); ix ++)
{
 for(vector<int>::size_type iy = 0; iy != pbonds.size(); iy += 2){
    if (ptetra[ix] == pbonds[iy])
    {
        pta.push_back(pbonds[iy+1]);
    }
}
for(vector<int>::size_type iy = 1; iy != (pbonds.size() + 1); iy += 2){
    if (ptetra[ix] == pbonds[iy])
    {
        pta.push_back(pbonds[iy-1]);
    }
}
}

}

/*
void print_to_file(ofstream myfile)
{
        // out put number of atoms
        myfile << "Set i /a1 * a"<<(*m).NumAtoms() + (*m1).NumAtoms()<<"/;\n";
        myfile << "Alias(i,j,k,l,m,n,o,p);\n";

        myfile << "* reactants\n";

    // output all atom types
        myfile << "Parameter rtypes(i) 'types of atoms'/";
        if (rtypes.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rtypes.size(); ix++)
            {
                if (ix==rtypes.size() - 1)
                {
                    myfile << "a"<< ix+1 <<" "<<rtypes[ix]<<"/;\n";
                } else {
                    myfile << "a"<< ix+1 <<" "<<rtypes[ix]<<", ";
                }   
            }
        } else {
            myfile << "/;\n";
        }

    // output all bond connection between two atoms for gams file
        myfile << "Set rbonds(i,j) 'bonds btween two atoms' /";

        if (rbonds.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rbonds.size(); ix += 2)
            {
                if (ix==rbonds.size() - 2)
                {
                    myfile << "a"<<rbonds[ix]<<"."<< "a"<<rbonds[ix+1]<<"/;\n";
                } else {
                    myfile << "a"<<rbonds[ix]<<"."<< "a"<<rbonds[ix+1]<<", ";
                }   
            }
        } else {
            myfile << "/;\n";
        }

// output tetra atoms
        myfile << "Set adjr /t1 * t4/;\n";
        myfile << "Alias(adjr,adjp);\n";
        myfile << "Set rtetra(i) 'index of tetra atom' /";
        if (rtetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rtetra.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< rtetra[ix] <<" /;\n";
                } else {
                    myfile << "a"<< rtetra[ix] <<" , ";
                }   
            }
        } else {
            myfile << "/;\n";
        }


// output four neighoring atoms of each tetra atoms
        myfile << "Parameter rta(i,adjr) 'neighoring atoms of each tetra atoms'/";
        if (rtetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rtetra.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< rtetra[ix] <<".t1 "<< rta[4*ix]<<" , "
                    << "a"<< rtetra[ix] <<".t2 "<< rta[4*ix+1]<<" , "
                    << "a"<< rtetra[ix] <<".t3 "<< rta[4*ix+2]<<" , "
                    << "a"<< rtetra[ix] <<".t4 "<< rta[4*ix+3]<<" /;\n";
                } else {
                    myfile << "a"<< rtetra[ix] <<".t1 "<< rta[4*ix]<<" , "
                    << "a"<< rtetra[ix] <<".t2 "<< rta[4*ix+1]<<" , "
                    << "a"<< rtetra[ix] <<".t3 "<< rta[4*ix+2]<<" , "
                    << "a"<< rtetra[ix] <<".t4 "<< rta[4*ix+3]<<" , ";
                }   
            }
        } else {
            myfile << "/;\n";
        }

// output types of tetra atoms for gams file
        myfile << "Set rtt(i) 'types of tetra atoms' /";
        if (rtetra.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rtetra.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< rtetra[ix] <<" 6/;\n";
                } else {
                    myfile << "a"<< rtetra[ix] <<" 6, ";
                }   
            }
        } else {
            myfile << "/;\n";
        }

// output double bonds for gams file
        myfile << "Set rdb(i,j) 'double bonds between two atom index'/";
        if (rdb.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rdb.size(); ix += 2)
            {
                if (ix==rdb.size() - 2)
                {
                    myfile << "a"<<rdb[ix]<<"."<< "a"<<rdb[ix+1]<<"/;\n";
                }
                myfile << "a"<<rdb[ix]<<"."<< "a"<<rdb[ix+1]<<", ";   
            }
        } else {
            myfile << "/;\n";
        }


// output four neighoring atoms of double bonds atoms
        myfile << "Parameter rda(i,j,adjr) 'neighoring atoms of double bonds'/";
        if (rdb.size()>1)
        {
            for(vector<int>::size_type ix = 0; ix != rdb.size(); ix++)
            {
                if (ix==rtetra.size() - 1)
                {
                    myfile << "a"<< rdb[ix] <<".t1 "<< rta[4*ix]<<" , "
                    << "a"<< rdb[ix] <<".t2 "<< rta[4*ix+1]<<" /;\n";
                } else {
                    myfile << "a"<< rdb[ix] <<".t1 "<< rta[4*ix]<<" , "
                    << "a"<< rdb[ix] <<".t2 "<< rta[4*ix+1]<<" , ";
                }   
            }
        } else {
            myfile << "/;\n";
        }
}

*/


        /*
        // get information of atoms
        // atom index start with 1 in OpenBabel
        int num_atoms = (*m).NumAtoms();
        std::vector<int> rtypes;
        for (int index = 1; index <= num_atoms; index++)
        {
            OBAtom atom = *((*m).GetAtom(index));
            std::cout << "Type of atom "<< *(atom.GetType()) << atom.GetAtomicMass() << std::endl;
            rtypes.push_back(int(atom.GetAtomicNum()));

        }        
        
        // get information of bonds
        // bond index start with 0 in OpenBabel
        int num_bonds = (*m).NumBonds();
        std::cout << "Number of bonds "<< num_bonds << std::endl;
        std::vector<int> rbonds;
        std::vector<int> rdb;
        
        for (int index = 0; index < num_bonds; index++)
        {
            OBBond bond = *((*m).GetBond(index));
            std::cout << "bond atom"<< (*(bond.GetBeginAtom())).GetIdx() 
            << *((*(bond.GetBeginAtom())).GetType()) << " bond atom"
            << (*(bond.GetEndAtom())).GetIdx() << *((*(bond.GetEndAtom())).GetType())
            << std::endl;

            // bond connection between two atoms
            rbonds.push_back((*(bond.GetBeginAtom())).GetIdx());
            rbonds.push_back((*(bond.GetEndAtom())).GetIdx());


            // double bonds
            if (bond.IsDouble()) 
            {
                rdb.push_back((*(bond.GetBeginAtom())).GetIdx());
                rdb.push_back((*(bond.GetEndAtom())).GetIdx());
            }
        }

        // determine which is tetra atom
        std::vector<int> rtetra;
        //std::vector<int> rtt;
        for (int i = 0; i <= num_atoms ; ++i)
        {
            int tetra_count = 0;
            if (rtypes[i] == 6) 
            {
                for(vector<int>::size_type ix = 0; ix != rbonds.size(); ix ++){
                    if (rbonds[ix] == (i + 1)) 
                    {
                        tetra_count += 1;
                    }
                }
                if (tetra_count == 4)
                {
                   rtetra.push_back(i+1);
               } 
           }
       }

     // find four neighoring atoms to double bonds
       std::vector<int> rda;
       for (vector<int>::size_type ix = 0; ix != rdb.size(); ix ++)
       {
         for(vector<int>::size_type iy = 0; iy != rbonds.size(); iy += 2){
            if (rdb[ix] == rbonds[iy] && rdb[ix+1] != rbonds[iy])
            {
                rda.push_back(rbonds[iy+1]);
            }
        }
        for(vector<int>::size_type iy = 1; iy != (rbonds.size() + 1); iy += 2){
            if (rdb[ix] == rbonds[iy] && rdb[ix+1] != rbonds[iy])
            {
                rda.push_back(rbonds[iy-1]);
            }
        }
    }

    // find four neighoring atoms to tetra atom
    std::vector<int> rta;
    for (vector<int>::size_type ix = 0; ix != rtetra.size(); ix ++)
    {
     for(vector<int>::size_type iy = 0; iy != rbonds.size(); iy += 2){
        if (rtetra[ix] == rbonds[iy])
        {
            rta.push_back(rbonds[iy+1]);
        }
    }
    for(vector<int>::size_type iy = 1; iy != (rbonds.size() + 1); iy += 2){
        if (rtetra[ix] == rbonds[iy])
        {
            rta.push_back(rbonds[iy-1]);
        }
    }
}

*/

    /*if(argc<3)
    {
     cout << "Usage: ProgrameName InputFileName OutputFileName\n";
     return 1;
    }
    ifstream ifs(argv[1]);
    if(!ifs)
    {
     cout << "Cannot open input file\n";
     return 1;
    }
    ofstream ofs(argv[2]);
    if(!ofs)
    {
     cout << "Cannot open output file\n";
     return 1;
    }
    OpenBabel::OBConversion conv;
    conv.SetInFormat("sdf");
    OpenBabel::OBConversion conv(&ifs, &ofs);
    conv(&ifs, &ofs);
    if(!conv.SetInAndOutFormats("CML","MOL"))
    {
     cout << "Formats not available\n";
     return 1;
    }
    int n = conv.Convert();
    cout << n << " molecules converted\n";
    return 0;*/
