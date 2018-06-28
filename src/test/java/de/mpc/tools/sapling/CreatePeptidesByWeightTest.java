package de.mpc.tools.sapling;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.junit.Test;

import de.mpc.tools.sapling.CreatePeptidesByWeight;
import uk.ac.ebi.pride.utilities.mol.AminoAcidSequence;

public class CreatePeptidesByWeightTest {

    @Test
    public void createRandomPeptide() {
        HashMap<Character, Double> fixedMods = new HashMap<Character, Double>();
        fixedMods.put('C', 57.021464);

        HashMap<Character, List<Double>> variableMods = new HashMap<Character, List<Double>>();
        variableMods.put('M', new ArrayList<Double>());
        variableMods.get('M').add(15.994915);

        CreatePeptidesByWeight pepCreator = new CreatePeptidesByWeight(fixedMods, variableMods, 12345L);

        for (int i = 0; i < 100; i++) {
            // TODO: useful testing...
            AminoAcidSequence sequence = pepCreator.createRandomPeptide(2, 2, 445.770996, 5.0);

            System.out.println("pep: " + sequence + ", weight: " + sequence.getMonoMass());
        }
    }
}
