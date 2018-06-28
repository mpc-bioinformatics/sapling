package de.mpc.tools.sapling;

import static uk.ac.ebi.pride.utilities.mol.NeutralLoss.WATER_LOSS;
import static uk.ac.ebi.pride.utilities.mol.NuclearParticle.PROTON;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import uk.ac.ebi.pride.utilities.mol.AminoAcid;
import uk.ac.ebi.pride.utilities.mol.AminoAcidSequence;
import uk.ac.ebi.pride.utilities.mol.MoleculeUtilities;

public class CreatePeptidesByWeight {

    /** the available aminoacids (no unknown, multiple or selenocysteine) */
    public static final List<AminoAcid> aminoAcidList = Collections.unmodifiableList(Arrays
            .asList(new AminoAcid[] {
                    AminoAcid.A, AminoAcid.R, AminoAcid.N, AminoAcid.D, AminoAcid.C, AminoAcid.E,
                    AminoAcid.Q, AminoAcid.G, AminoAcid.H, AminoAcid.I, AminoAcid.L, AminoAcid.K,
                    AminoAcid.M, AminoAcid.F, AminoAcid.P, AminoAcid.S, AminoAcid.T, AminoAcid.V,
                    AminoAcid.W, AminoAcid.Y }));

    /** the number of available AAs */
    public static final int nrAAs = aminoAcidList.size();

    /** random number generator */
    private Random rnd;

    /** fixed modifications, mapping from the single letter AA to weight */
    private HashMap<Character, Double> fixedMods;

    /**
     * variable modifications, mapping from the single letter AA to possible
     * weights
     */
    private HashMap<Character, List<Double>> variableMods;

    private static double[][] aaDiffs;

    public CreatePeptidesByWeight(HashMap<Character, Double> fixedMods, HashMap<Character, List<Double>> variableMods) {
        this.fixedMods = (fixedMods == null) ? new HashMap<Character, Double>() : fixedMods;
        this.variableMods = (variableMods == null) ? new HashMap<Character, List<Double>>() : variableMods;

        initializeDifferences();
    }

    public CreatePeptidesByWeight(HashMap<Character, Double> fixedMods, HashMap<Character, List<Double>> variableMods,
            long seed) {
        this(fixedMods, variableMods);
        rnd = new Random(seed);
    }

    /**
     * Calculate the differences between amino acids, taking fixed modifications
     * into account.
     */
    private void initializeDifferences() {
        aaDiffs = new double[nrAAs][nrAAs];

        for (int i = 0; i < nrAAs; i++) {
            double massI = aminoAcidList.get(i).getMonoMass();
            if (fixedMods.containsKey(aminoAcidList.get(i).getOneLetterCode())) {
                massI += fixedMods.get(aminoAcidList.get(i).getOneLetterCode());
            }

            for (int j = 0; j < nrAAs; j++) {
                double massJ = aminoAcidList.get(j).getMonoMass();
                if (fixedMods.containsKey(aminoAcidList.get(j).getOneLetterCode())) {
                    massJ += fixedMods.get(aminoAcidList.get(j).getOneLetterCode());
                }

                aaDiffs[i][j] = massI - massJ;
            }
        }
    }

    /**
     * Tries to create a random AA sequence, that fits into the given charge,
     * m/z and deltaPPM values. If none could be after several retries, null
     * will be returned.
     *
     * @param minCharge
     * @param maxCharge
     * @param precursorMz
     * @param maxDeltaPPM
     * @return either an appropriate AA sequence or null, if none could be found
     */
    public AminoAcidSequence createRandomPeptide(int minCharge, int maxCharge, double precursorMz, double maxDeltaPPM) {
        int tries = 0;

        AminoAcidSequence sequence = null;
        while ((tries < 1000) && (sequence == null)) {
            sequence = createRandomPeptide(rnd.nextInt(maxCharge - minCharge + 1) + minCharge, precursorMz,
                    maxDeltaPPM);
        }

        return sequence;
    }

    /**
     * Tries to create a random AA sequence, that fits into the given charge,
     * m/z and deltaPPM values.
     *
     * @param precursorCharge
     * @param precursorMz
     * @param maxDeltaPPM
     * @return either an appropriate AA sequence or null, if none could be found
     */
    private AminoAcidSequence createRandomPeptide(int precursorCharge, double precursorMz, double maxDeltaPPM) {
        // calculate the precursor mass from the precursorMz and charge
        double precursorMass = Math.abs(precursorMz * precursorCharge - precursorCharge * PROTON.getMonoMass());

        // the variable modifications
        ArrayList<Double> vMods = new ArrayList<Double>();

        double deltaPPM;

        AminoAcidSequence sequence = new AminoAcidSequence();
        boolean inWeight = true;
        while (inWeight) {
            AminoAcid aa = aminoAcidList.get(rnd.nextInt(nrAAs));

            AminoAcidSequence newSequence = new AminoAcidSequence();
            newSequence.addAminoAcids(sequence.getAminoAcids());
            newSequence.addAminoAcid(aa);

            if ((variableMods.containsKey(aa.getOneLetterCode())) && rnd.nextBoolean()) {
                List<Double> mods = variableMods.get(aa.getOneLetterCode());
                vMods.add(mods.get(rnd.nextInt(mods.size())));
            } else {
                vMods.add(0.0);
            }

            deltaPPM = calculatePPM(newSequence, vMods, precursorCharge, precursorMass);

            if (deltaPPM > -maxDeltaPPM) {
                sequence = newSequence;
            }

            // check, whether adding another AA is possible for the allowed MZ
            // and PPM (checking against glycine)
            newSequence = new AminoAcidSequence();
            newSequence.addAminoAcids(sequence.getAminoAcids());
            newSequence.addAminoAcid(AminoAcid.G);

            deltaPPM = calculatePPM(newSequence, vMods, precursorCharge, precursorMass);

            if (deltaPPM < -maxDeltaPPM) {
                inWeight = false;
            }
        }

        // adding another AA makes peptide too heavy -> exchanging some to fit
        // in
        deltaPPM = calculatePPM(sequence, vMods, precursorCharge, precursorMass);

        int tries = 0;
        while (((deltaPPM > maxDeltaPPM) || (deltaPPM < -maxDeltaPPM)) && (tries < 100)) {
            double theoreticalMass = MoleculeUtilities.calculateTheoreticalMass(sequence.getOneLetterCodeString(),
                    createMassesArrayWithFixed(sequence, vMods));
            double diff = precursorMass - theoreticalMass;

            int bestI = 0;
            int bestJ = 0;
            double bestDiff = Double.POSITIVE_INFINITY;

            // looking for best fit
            for (int i = 0; i < nrAAs; i++) {
                for (int j = 0; j < nrAAs; j++) {
                    double dDiff = diff - aaDiffs[i][j];

                    if (Math.abs(dDiff) < bestDiff) {
                        bestI = i;
                        bestJ = j;
                        bestDiff = Math.abs(dDiff);
                    }

                    // TODO: add variable modifications here at the moment,they
                    // are ignored, which might create a suboptimal swapping
                }
            }

            if ((bestJ != bestI) && (sequence.getAminoAcids().contains(aminoAcidList.get(bestJ)))) {
                // exchange the first occurrence of the amino acid
                String seq = sequence.getOneLetterCodeString();
                int pos = seq.indexOf(aminoAcidList.get(bestJ).toString());

                // System.out.println("\t" + sequence + "\t swapping " +
                // aminoAcidList.get(bestJ) + " -> " +
                // aminoAcidList.get(bestI));

                sequence.remove(pos);
                sequence.addAminoAcid(aminoAcidList.get(bestI));

                vMods.remove(pos);
                vMods.add(0.0);
            } else {
                int pos = rnd.nextInt(sequence.getLength());

                // System.out.println("\t" + sequence + "\t random remove of " +
                // sequence.getAminoAcid(pos) + " at " + pos);

                sequence.remove(pos);
                sequence.addAminoAcid(aminoAcidList.get(rnd.nextInt(nrAAs)));

                vMods.remove(pos);
                vMods.add(0.0);
            }

            deltaPPM = calculatePPM(sequence, vMods, precursorCharge, precursorMass);

            tries++;
        }

        // -----------------------------
        /*
         * StringBuilder sb = new StringBuilder(sequence.getLength()); for (int
         * i=0; i<sequence.getLength(); i++) {
         *
         * sb.append( (vMods.get(i) != 0.0) ?
         * Character.toLowerCase(sequence.getAminoAcid(i).getOneLetterCode()) :
         * sequence.getAminoAcid(i).getOneLetterCode() ); }
         *
         * System.out.println(sb.toString() + "\t c: " + precursorCharge +
         * "\t ppm: " + deltaPPM + "\t tries: " + tries);
         */
        // -----------------------------

        if (Math.abs(deltaPPM) < maxDeltaPPM) {
            return sequence;
        } else {
            return null;
        }
    }

    /**
     *
     * @param masses
     * @return
     */
    private double[] createMassesArrayWithFixed(AminoAcidSequence sequence, List<Double> vMods) {
        ArrayList<Double> massShifts = new ArrayList<Double>();
        // always has water loss due to the enzyme
        massShifts.add(WATER_LOSS.getMonoMass());

        // add the masses for fixed modifications
        for (Character aa : sequence.getOneLetterCodeString().toCharArray()) {
            if (fixedMods.containsKey(aa)) {
                massShifts.add(fixedMods.get(aa));
            }
        }

        // add the masses of the variable modifications
        for (Double mass : vMods) {
            if (!mass.equals(0.0)) {
                massShifts.add(mass);
            }
        }

        double[] massesArray = new double[massShifts.size()];
        for (int i = 0; i < massShifts.size(); i++) {
            massesArray[i] = massShifts.get(i);
        }

        return massesArray;
    }

    /**
     * Calculates the PPM shift of the given sequence.
     *
     * @param sequence
     * @param vMods
     *            the used, variable modifications
     * @param precursorCharge
     * @param precursorMass
     * @return
     */
    private double calculatePPM(AminoAcidSequence sequence, List<Double> vMods, int precursorCharge,
            double precursorMass) {
        double theoreticalMass = MoleculeUtilities.calculateTheoreticalMass(sequence.getOneLetterCodeString(),
                createMassesArrayWithFixed(sequence, vMods));
        double theoreticalMz = theoreticalMass / precursorCharge + PROTON.getMonoMass();
        double deltaMZ = (precursorMass - theoreticalMass) / precursorCharge;

        return deltaMZ / theoreticalMz * 1000000.0;
    }
}
