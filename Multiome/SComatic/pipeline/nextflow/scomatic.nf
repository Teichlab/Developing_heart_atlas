nextflow.enable.dsl=2

// These are the expected command line arguments
// GEX/ATAC configuration to be delivered via a -params-file
params.mappings = null
params.celltypes = null
// We expect the mappings to be on irods by default
params.location = "irods"
// The genotyping pipeline will specify an scomatic mutations folder
params.mutations = null

// Download a given sample's BAM from iRODS
// Then either retrieve the BAI or make one via indexing
// The maxForks of 10 was set after asking jc18 about best iRODS practices
process irods {
    maxForks 10
    label "normal4core"
    input:
        tuple val(sample), val(irods), val(bam), val(donor)
    output:
        tuple val(sample), val(donor), path("${sample}.bam"), path("${sample}.bam.bai")
    script:
        """
        iget -K ${irods}/${bam} ${sample}.bam
        if [[ `ils ${irods} | grep "${bam}.bai" | wc -l` == 1 ]]
        then
            iget -K ${irods}/${bam}.bai ${sample}.bam.bai
        else
            samtools index -@ ${task.cpus} ${sample}.bam
        fi
        """
}

// The equivalent of an irods download, but for a local copy of mappings
// Symlink the BAM/BAI appropriately so they're named the right thing for downstream
process local {
    label "normal4core"
    input:
        tuple val(sample), path(local), val(bam), val(donor)
    output:
        tuple val(sample), val(donor), path("${sample}.bam"), path("${sample}.bam.bai")
    script:
        """
        ln -s ${local}/${bam} ${sample}.bam
        if [ -f "${local}/${bam}.bai" ]
        then
            ln -s ${local}/${bam}.bai ${sample}.bam.bai
        else
            samtools index -@ ${task.cpus} ${sample}.bam
        fi
        """
}

// Extract the cell type definitions for the specified sample
// Returns the sample as part of the tuple for easy joining with sample BAM/BAI lists
// Finds the sample ID at the start of the line as per SAMPLE_BARCODE-1 specification
process getSampleCelltypes {
    label "normal"
    input:
        tuple val(sample), path(celltypes)
    output:
        tuple val(sample), path("${sample}-celltypes.tsv")
    script:
        """
        head -n 1 ${celltypes} > ${sample}-celltypes.tsv
        grep "^${sample}" ${celltypes} >> ${sample}-celltypes.tsv
        """
}

// Step one of scomatic - split the BAM into per cell type BAMs
// However, in the interest of efficiency, we do this per sample
// So we have to merge them on a donor level across samples later
// GEX/ATAC alternate parameterisation (NH/nM vs MQ) prepared in modargs
// The output for this is a little weird, as the process makes a "list of lists"
// Meanwhile I need to create a nice clean tuple for subsequent alternate grouping
// As such, this emits BAM/BAI "list of lists", which have all the necessary info
// (sample.donor.celltype) embedded in the file names
// The typical scomatic convention is sample.celltype, so I have to rename files
process splitBams {
    label "long"
    input:
        tuple val(sample), val(donor), path(bam), path(bai), path(celltypes), val(modargs)
    output:
        path("output/*.bam"), emit: "bam"
        path("output/*.bam.bai"), emit: "bai"
    script:
        """
        mkdir -p output
        python3 ${params.scomatic}/SplitBam/SplitBamCellTypes.py --bam ${bam} \
            --meta ${celltypes} \
            --id ${sample} \
            --n_trim 5 \
            ${modargs} \
            --outdir output
        rename "s/${sample}/${sample}.${donor}/g" output/*
        """
}

// Merge the per-sample BAMs for each cell type into one per-donor file
// Pass the donor and cell type along for a while for downstream process use
process mergeCelltypeBams {
    label "normal4core"
    input:
        tuple val(donor), val(celltype), path(bams), path(bais)
    output:
        tuple val(donor), val(celltype), path("${donor}.${celltype}.bam")
    script:
        """
        samtools merge -@ {task.cpus} ${donor}.${celltype}.bam ${bams}
        """
}

// Index the donor-celltype BAMs, returning both the BAM and BAI
process indexCelltypeBams {
    label "normal4core"
    input:
        tuple val(donor), val(celltype), path(bam)
    output:
        tuple val(donor), val(celltype), path(bam), path("*.bam.bai")
    script:
        """
        samtools index -@ {task.cpus} ${bam}
        """
}

// Step two of scomatic - convert the donor-celltype BAMs to scomatic's TSVs
// Up the bin parameter to decrease the number of temporary files made
// Also the LSF config for this uses scratch to not pollute the drive with those files
// (Even if they are just temporary and get removed at the end of the process)
// There can sometimes be no output file, so specify it as optional
// Can't pass out donor info because tuples and optional don't get along
process bamToTsv {
    cpus 16
    label "long16core"
    input:
        tuple val(donor), val(celltype), path(bam), path(bai), path(fasta), path(fai), val(mq)
    output:
        path("output/*.tsv", optional: true)
    script:
        """
        mkdir -p temp
        mkdir -p output
        python3 ${params.scomatic}/BaseCellCounter/BaseCellCounter.py --bam ${bam} \
            --ref ${fasta} \
            --chrom all \
            --out_folder output \
            --min_bq 30 \
            --min_mq ${mq} \
            --tmp_dir temp \
            --bin 1000000 \
            --nprocs ${task.cpus}
        rm -rf temp
        """
}

// Step three of scomatic - merge the celltype TSVs into a single mega TSV
// Important to stage the input in a directory here so scomatic sees it like it wants to
process mergeTsvs {
    label "week"
    input:
        tuple val(donor), path(tsvs, stageAs: "input/*")
    output:
        tuple val(donor), path("${donor}.BaseCellCounts.AllCellTypes.tsv")
    script:
        """
        python3 ${params.scomatic}/MergeCounts/MergeBaseCellCounts.py --tsv_folder input \
            --outfile ${donor}.BaseCellCounts.AllCellTypes.tsv
        """
}

// Step 4.1 of scomatic - perform the initial mutation calling
// Set max_cell_types to a stratospheric value to disable that filter
// As we are interested in germline mutations spanning multiple populations
process callMutations {
    label "week"
    input:
        tuple val(donor), path(tsv), path(fasta), path(fai)
    output:
        tuple val(donor), path("${donor}.calling.step1.tsv")
    script:
        """
        python3 ${params.scomatic}/BaseCellCalling/BaseCellCalling.step1.py \
            --infile ${tsv} \
            --outfile ${donor} \
            --ref ${fasta} \
            --max_cell_types 1000000
        """
}

// Step 4.2 of scomatic - perform initial mutation filtering
// Need to do GEX and ATAC separately due to the editing file
// This needs a bunch of memory apparently for reasons that I can't identify from the code
process filterMutationsGex {
    label "long10gb"
    publishDir "${donor}-${params.modality}", mode:"copy"
    input:
        tuple val(donor), path(tsv), path(pons), path(editing)
    output:
        tuple val(donor), path("${donor}.calling.step2.tsv")
    script:
        """
        python3 ${params.scomatic}/BaseCellCalling/BaseCellCalling.step2.py \
            --infile ${tsv} \
            --outfile ${donor} \
            --editing ${editing} \
            --pon ${pons}
        """
}

process filterMutationsAtac {
    label "long10gb"
    publishDir "${donor}-${params.modality}", mode:"copy"
    input:
        tuple val(donor), path(tsv), path(pons)
    output:
        tuple val(donor), path("${donor}.calling.step2.tsv")
    script:
        """
        python3 ${params.scomatic}/BaseCellCalling/BaseCellCalling.step2.py \
            --infile ${tsv} \
            --outfile ${donor} \
            --pon ${pons}
        """
}

// Intersect the filtered mutations with a BED region of interest
process intersectBed {
    label "normal"
    publishDir "${donor}-${params.modality}", mode:"copy"
    input:
        tuple val(donor), path(tsv), path(bed)
    output:
        tuple val(donor), path("${donor}.calling.step2.intersect.tsv")
    script:
        """
        bedtools intersect -header -a ${tsv} -b ${bed} > ${donor}.calling.step2.intersect.tsv
        """
}

// Final mutation filter, just keep the PASS ones
// Standard $ escaping applies
process passMutations {
    label "normal"
    publishDir "${donor}-${params.modality}", mode:"copy"
    input:
        tuple val(donor), path(tsv)
    output:
        tuple val(donor), path("${donor}.calling.step2.pass.tsv")
    script:
        """
        awk '\$1 ~ /^#/ || \$6 == "PASS"' ${tsv} > ${donor}.calling.step2.pass.tsv
        """
}

// Get callable sites on a per cell type level
// This makes a couple files, both end in .report.tsv
process callableSitesCellType {
    label "long"
    publishDir "${donor}-${params.modality}", mode:"copy"
    input:
        tuple val(donor), path(tsv)
    output:
        tuple val(donor), path("*.report.tsv")
    script:
        """
        python3 ${params.scomatic}/GetCallableSites/GetAllCallableSites.py \
            --infile ${tsv} \
            --outfile ${donor}
        """
}

// Get callable sites on a cell level
// Mirror argument values based on earlier bamToTsv step for consistency
process callableSitesCell {
    label "long16core"
    publishDir "${donor}-${params.modality}/cell_callable_sites", mode:"copy"
    input:
        tuple val(donor), val(celltype), path(bam), path(bai), path(fasta), path(fai), val(mq), path(tsv)
    output:
        tuple val(donor), val(celltype), path("*.SitesPerCell.tsv")
    script:
        """
        mkdir -p temp
        python3 ${params.scomatic}/SitesPerCell/SitesPerCell.py --bam ${bam} \
            --infile ${tsv} \
            --ref ${fasta} \
            --min_bq 30 \
            --min_mq ${mq} \
            --tmp_dir temp \
            --bin 1000000 \
            --nprocs ${task.cpus}
        rm -rf temp
        """
}

// Post-processing - convert the donor-celltype BAMs to per cell genotypes
// Up the bin parameter to decrease the number of temporary files made
// Also the LSF config for this uses scratch to not pollute the drive with those files
// (Even if they are just temporary and get removed at the end of the process)
process bamToGenotype {
    label "long16core10gb"
    publishDir "${donor}-${params.modality}-genotypes", mode:"copy"
    input:
        tuple val(donor), val(celltype), path(bam), path(bai), path(fasta), path(fai), val(mq), path(allcelltypes), path(mutations)
    output:
        path("${donor}.${celltype}.single_cell_genotype.tsv", optional: true)
    script:
        """
        mkdir -p temp
        python3 ${params.scomatic}/SingleCellGenotype/SingleCellGenotype.py --bam ${bam} \
            --ref ${fasta} \
            --infile ${mutations}/${donor}-${params.modality}/${donor}.calling.step2.intersect.tsv \
            --meta ${allcelltypes} \
            --outfile ${donor}.${celltype}.single_cell_genotype.tsv \
            --min_bq 30 \
            --min_mq ${mq} \
            --tmp_dir temp \
            --bin 1000000 \
            --nprocs ${task.cpus}
        rm -rf temp
        """
}

// Perform the various data downloading and demultiplexing early in scomatic
// Yields demultiplexed cell type BAMs on a donor level
// Along with the master cell type file and the reference genome + index
// And the upcoming downstream process --min_mq (255 GEX, 30 ATAC)
// As both downstream uses want that to call the process appropriately
// The cell type BAMs are a good input for both normal scomatic
// And the optional secondary pipeline of single cell genotype calling
workflow STEP1 {
    main:
        // Make sure all the input files actually exist
        if (params.mappings == null) {
            error "Please provide a mappings CSV file via --mappings"
        }
        if (params.celltypes == null) {
            error "Please provide a cell types TSV file via --celltypes"
        }
        // Load the mapping as a CSV, skipping the one line of header for ease of use
        // Split the BAM file into its parent directory and the actual file name
        // As the parent directory is used both by irods and local for finding the BAI
        // It's fine to use file() here for the file name operations
        // Nextflow doesn't natively check the paths exist, allowing for irods splitting
        mappings = Channel
            .fromPath(params.mappings, checkIfExists: true)
            .splitCsv(skip:1)
            .map({it -> [it[0], file(it[1]).getParent(), file(it[1]).getName(), it[2]]})
        allCelltypes = Channel.fromPath(params.celltypes, checkIfExists: true)
        fasta = Channel.fromPath(params.genome, checkIfExists: true)
        fai = Channel.fromPath(params.genome+".fai", checkIfExists: true)
        // Check the modality
        if (!(params.modality in ["GEX","ATAC"])) {
            error "Unknown modality, must be GEX or ATAC"
        }
        // Sanitise the mappings to just samples that show up in the cell types file
        // To do this, we need a list of sample IDs from the cell type file
        // Need to use splitCsv as regular splitText doesn't have a skip option
        // Skip the header line, then keep only the first column for the file
        // And split that string on _ and remove the last element from it
        // Then glue it back together with _ as the delimiter
        // The equivalent of "_".join(VARIABLE.split("_")[:-1]) in python
        celltypeSamples = allCelltypes
            .splitCsv(sep:"\t", skip:1)
            .map({it -> it[0].split("_")[0..-2].join("_")})
            .unique()
        // Doing a join here keeps the intersection of the sample IDs
        // That show up in both the mappings info and our parsed celltype sample list
        mappings = mappings.join(celltypeSamples)
        if (params.location == "irods") {
            // Download all the mappings from irods - yields a sample/donor/BAM/BAI tuple
            sampleBams = irods(mappings)
        }
        else if (params.location == "local") {
            // The mappings are available locally, make the same tuple
            // Symlink up the files so they're named like the downloads would have been
            sampleBams = local(mappings)
        }
        else {
            error "Unknown location, must be irods or local"
        }
        // Extract the samples as the first element of the tuple
        // Then combine it with the cell types file
        samples = sampleBams
            .map({it -> it[0]})
            .combine(allCelltypes)
        // Can now search the cell types file for each of the sample IDs
        sampleCelltypes = getSampleCelltypes(samples)
        // Join the cell type file with the BAM/BAI list from earlier
        sampleBamsCelltypes = sampleBams.join(sampleCelltypes)
        // Prepare modality-dependent arguments as an extra argument to the process
        if (params.modality == "GEX") {
            sampleBamsCelltypes = sampleBamsCelltypes.combine(Channel.of("--max_nM 5 --max_NH 1"))
        } else if (params.modality == "ATAC") {
            sampleBamsCelltypes = sampleBamsCelltypes.combine(Channel.of("--min_MQ 30"))
        }
        // Step one of scomatic - split the BAMs to cell types on a per sample basis
        sampleFiles = splitBams(sampleBamsCelltypes)
        // This outputs a "list of lists" of BAMs/BAIs, constructed for each sample
        // Flatten the files into a single list, and then index them as 
        // sample.donor.celltype
        // (the file names are sample.donor.celltype.bam[.bai])
        // Note that the dot needs to be escaped for the parsing to work, here and later
        bam = sampleFiles.bam
            .flatten()
            .map({file -> [file.getName().split('\\.bam')[0], file]})
        bai = sampleFiles.bai
            .flatten()
            .map({file -> [file.getName().split('\\.bam')[0], file]})
        // We can now safely join the two file lists on the index to combine correctly
        // Replace the index with two entries, one for donor, one for cell type
        // And then group up all the files from a given donor and cell type for merging
        sampleCelltypeBams = bam
            .join(bai)
            .map({name, file, ind -> [name.split('\\.')[1], name.split('\\.')[2], file, ind]})
            .groupTuple(by: [0,1])
        unindexedCelltypeBams = mergeCelltypeBams(sampleCelltypeBams)
        // Index the cell type BAMs, and then add the genome to the resulting tuple
        // All downstream uses of the BAMs want the genome present, so store that too
        indexedCelltypeBams = indexCelltypeBams(unindexedCelltypeBams)
            .combine(fasta)
            .combine(fai)
        // Let's also store the upcoming min_mq value as used by all processes
        // 255 for GEX, 30 for ATAC
        if (params.modality == "GEX") {
            indexedCelltypeBams = indexedCelltypeBams.combine(Channel.of("255"))
        } else if (params.modality == "ATAC") {
            indexedCelltypeBams = indexedCelltypeBams.combine(Channel.of("30"))
        }
    emit:
        indexedCelltypeBams = indexedCelltypeBams
        allCelltypes = allCelltypes
        fasta = fasta
        fai = fai
}

// Turn an existing set of scomatic mutations, plus freshly re-generated cell type BAMs
// Into single cell level mutation calls
workflow genotypes {
    // A folder with prior scomatic mutations output
    // Needs DONOR-MODALITY subfolders present
    if (params.mutations == null) {
        error "Please provide a path to prior scomatic pipeline mutation output via --mutations"
    }
    mutations = Channel.fromPath(params.mutations, checkIfExists: true)
    // Run the first step of scomatic, getting the cell type BAMs
    STEP1()
    // Pull out the various outputs that we use
    indexedCelltypeBams = STEP1.out.indexedCelltypeBams
    allCelltypes = STEP1.out.allCelltypes
    // We have the donor cell type genotypes we're after here
    // Just need to include the master cell types file and the master mutations folder
    indexedCelltypeBams = indexedCelltypeBams
        .combine(allCelltypes)
        .combine(mutations)
    // Get the single cell genotypes
    cellGenotypes = bamToGenotype(indexedCelltypeBams)
}

// The main scomatic workflow
workflow {
    // Some inputs that do not matter to step one, and as such are absent there
    // But we care about them here, so load them up and check them
    pons = Channel.fromPath(params.pons, checkIfExists: true)
    // This is a GEX-only file
    if (params.modality == "GEX") {
        editing = Channel.fromPath(params.editing, checkIfExists: true)
    }
    bed = Channel.fromPath(params.bed, checkIfExists: true)
    // Run the first step of scomatic, getting the cell type BAMs
    STEP1()
    // Pull out the various outputs that we use
    indexedCelltypeBams = STEP1.out.indexedCelltypeBams
    fasta = STEP1.out.fasta
    fai = STEP1.out.fai
    // Step two of scomatic - convert the BAMs to scomatic's TSV format for each cell type
    celltypeTsvs = bamToTsv(indexedCelltypeBams)
    // This only stores output where the thing actually ran, don't have to filter on that
    // Have to rebuild donor info from the file name as tuples and optional don't work
    // Then group them up on donor for collapsing to a per-donor TSV
    donorCelltypeTsvs = celltypeTsvs
        .map({file -> [file.getName().split('\\.')[0], file]})
        .groupTuple()
    // Step three of scomatic - merge the cell type TSVs into a single one per donor
    // We need to also pass the genome to the next scomatic step, so combine it in
    masterTsvs = mergeTsvs(donorCelltypeTsvs)
        .combine(fasta)
        .combine(fai)
    // Step 4.1 of scomatic - perform the initial mutation calling
    unfilteredMutationsClean = callMutations(masterTsvs)
    // Combine in the PoNs for filtering, but keep a clean version for later
    unfilteredMutations = unfilteredMutationsClean.combine(pons)
    // Step 4.2 of scomatic - filter the called mutations based on editing/PoNs
    // This requires two separate processes for the modalities as ATAC has no editing
    // Regardless of process, combine in the BED for later
    if (params.modality == "GEX") {
        filteredMutations = filterMutationsGex(unfilteredMutations.combine(editing))
            .combine(bed)
    }
    else if (params.modality == "ATAC") {
        filteredMutations = filterMutationsAtac(unfilteredMutations)
            .combine(bed)
    }
    // Final filtering steps - intersect with a BED of interest
    intersectedMutations = intersectBed(filteredMutations)
    // Stringently keep just the PASS mutations
    passedMutations = passMutations(intersectedMutations)
    // Compute the number of callable sites
    // Start off with per cell type, can just feed it the step 4.1 output
    cellTypeSites = callableSitesCellType(unfilteredMutationsClean)
    // The cell level callable sites require the cell type BAMs from earlier
    // To which we can "non-unique join" the step 4.1 output
    // The "non-unique join" is accomplished via .combine(by: 0)
    cellSitesInput = indexedCelltypeBams.combine(unfilteredMutationsClean, by: 0)
    // And now we can get cell level callable sites
    cellSites = callableSitesCell(cellSitesInput)
}