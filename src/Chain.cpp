//
// Created by jbarnett8 on 2/10/18.
//

#include "Chain.h"

using namespace std;
using namespace PNAB;
using namespace OpenBabel;

Chain::Chain(BaseUnit bu, unsigned chain_length, string ff_type) {

    vector<int> headIds, tailIds;
    OBAtom *head, *tail;
    monomer = bu.getMol();
    unsigned numMonomerAtoms = monomer.NumAtoms(), neighborId;
    chain_length_ = chain_length;

    chain = OBMol();
    auto bu_linkers = bu.getBackboneLinkers();


    // Creating raw chain
    for (unsigned polymerI = 0; polymerI < chain_length; polymerI++) {
        chain+=monomer;
    }

    // Deleting appropriate atoms and creating bonds
    for (int polymerI = 0; polymerI < chain_length - 1; polymerI++) {
        headIds.push_back(chain.GetAtom(numMonomerAtoms *   polymerI   + bu_linkers[0])->GetId());
        tailIds.push_back(chain.GetAtom(numMonomerAtoms * (polymerI+1) + bu_linkers[1])->GetId());
    }
    for (int polymerI = 0; polymerI < chain_length - 1; polymerI++) {
        head = chain.GetAtomById(headIds.at(polymerI));
        tail = chain.GetAtomById(tailIds.at(polymerI));
        FOR_NBORS_OF_ATOM(nbr, tail) {
            if (nbr->IsHydrogen()) {
                deletedAtomsId.push_back(nbr->GetId());
                chain.DeleteAtom(chain.GetAtom(nbr->GetIdx()));
                break;
            }
        }
        FOR_NBORS_OF_ATOM(nbr, head) {
            if (nbr->IsHydrogen()) {
                deletedAtomsId.push_back(nbr->GetId());
                chain.DeleteAtom(chain.GetAtom(nbr->GetIdx()));
                break;
            }
        }
        // chain.DeleteHydrogens(tail);  chain.DeleteHydrogens(head);
        FOR_NBORS_OF_ATOM(nbr, tail) {
            neighborId = nbr->GetId();
            break;
        }
        deletedAtomsId.push_back(tail->GetId());
        chain.DeleteAtom(tail);
        newBondIDs.push_back(head->GetId()); newBondIDs.push_back(neighborId);
        chain.AddBond(head->GetIdx(), chain.GetAtomById(neighborId)->GetIdx(), 1);
    }

    std::sort (deletedAtomsId.begin(), deletedAtomsId.end());

    pFF = OBForceField::FindForceField("GAFF");
    pFF->Setup(chain);
    isKCAL = pFF->GetUnit().find("kcal") != string::npos;

    // Setting up constraints
    constraintsTot = new OBFFConstraints;
    constraintsTor = new OBFFConstraints;
    constraintsAng = new OBFFConstraints;
    constraintsBond = new OBFFConstraints;
    vector<int> ignored_atoms;
    int nbrIdx;

    // Torsional energy, must fix non-torsion atoms
    vector<int> torsionAtoms;
    for (int i = 0; i < newBondIDs.size(); i++) {
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtomById(newBondIDs.at(i))) {
            torsionAtoms.push_back(nbr->GetIdx());
        }
    }
    FOR_ATOMS_OF_MOL(a, chain) {
        if(std::find(torsionAtoms.begin(), torsionAtoms.end(), a->GetIdx()) == torsionAtoms.end()) {
            constraintsTor->AddIgnore(a->GetIdx());
            constraintsAng->AddIgnore(a->GetIdx());
            constraintsBond->AddIgnore(a->GetIdx());
        }
    }

    // Angle energy, must fix non-angle atoms
    for (int i = 1; i < newBondIDs.size(); i+=2) {
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtomById(newBondIDs.at(i))) {
            nbrIdx = nbr->GetIdx();
            if (nbrIdx != chain.GetAtomById(newBondIDs.at(i - 1))->GetIdx()) {
                constraintsAng->AddIgnore(nbrIdx);
                constraintsBond->AddIgnore(nbrIdx);
            }
        }
    }

    // Bond energy, must fix non-bond atoms
    for (int i = 0; i < newBondIDs.size(); i+=2) {
        FOR_NBORS_OF_ATOM(nbr, chain.GetAtomById(newBondIDs.at(i))) {
            nbrIdx = nbr->GetIdx();
            if (nbrIdx != chain.GetAtomById(newBondIDs.at(i + 1))->GetIdx()) {
                constraintsBond->AddIgnore(nbrIdx);
            }
        }
    }

    if (!(constraintsTot && constraintsTor && constraintsAng && constraintsBond)) {
        cerr << "Error initializing force field constraints for chain. Exiting..." << endl;
        exit(1);
    }

}

ConformerData Chain::generateConformerData(double *conf, HelicalParameters hp) {
    if (chain.NumAtoms() <= 0) {
        cout << "There are no atoms in the chain. Exiting..." << endl;
        exit(0);
    }

    int numMonomerAtoms = monomer.NumAtoms(), ai, offset, dai = 0, xyzI = -3;
    double *xyz;
    auto g_rot = hp.getGlobalRotationOBMatrix();
    auto g_trans = hp.getGlobalTranslationVec();

    xyz = new double[chain.NumAtoms() * 3];

    for (int chainIdx = 0; chainIdx < chain_length_; chainIdx++) {
        for (int ci = 0; ci < 3*numMonomerAtoms; ci+=3) {
            offset = chainIdx * numMonomerAtoms;
            ai = (ci + 1) / 3 + offset + 1;
            /* for some reason, OpenBabel starts id's indexation at 0 instead of 1 like *
            *  the normal indexation for the atom index. This is the source of the +1  */
            if (dai == -1 || ai != deletedAtomsId.at(dai) + 1) {
                xyzI+=3;
                vector3 v;
                v.Set(conf + ci);
                v += hp.getStepTranslationVec(chainIdx);
                v *= hp.getStepRotationOBMatrix(chainIdx);
                v.Get(xyz + xyzI);
            } else {
                if ((dai + 1) < deletedAtomsId.size())
                    dai++;
                else
                    dai = -1;
            }
        }
    }

    ConformerData data;
    data.coords = xyz;
    data.chain_coords_present = true;
    fillConformerEnergyData(xyz, data);
    return data;
}

void Chain::fillConformerEnergyData(double *xyz, PNAB::ConformerData &conf_data) {
    chain.SetCoordinates(xyz);

    auto n = chain_length_;

    // Get total energies and VDW
    pFF->Setup(chain, *constraintsTot);
    conf_data.total_energy = pFF->Energy() / n;
    conf_data.VDWE = pFF->E_VDW() / n;
    conf_data.totTorsionE = pFF->E_Torsion() / n;

    if (n - 1 <= 0)
        n = 2;

    // Get torsional energies
    pFF->SetConstraints(*constraintsTor);
    conf_data.torsionE = pFF->E_Torsion() / (n - 1);

    // Get angle energy
    pFF->SetConstraints(*constraintsAng);
    conf_data.angleE = pFF->E_Angle() / (n - 1);

    // Get bond energy
    pFF->SetConstraints(*constraintsBond);
    conf_data.bondE = pFF->E_Bond() / (n - 1);

    if (!isKCAL) {
        conf_data.total_energy *= KJ_TO_KCAL;
        conf_data.VDWE         *= KJ_TO_KCAL;
        conf_data.totTorsionE  *= KJ_TO_KCAL;
        conf_data.torsionE     *= KJ_TO_KCAL;
        conf_data.angleE       *= KJ_TO_KCAL;
        conf_data.bondE        *= KJ_TO_KCAL;
    }
}
