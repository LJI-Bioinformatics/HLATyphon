import sys
import operator
import re
import os
import itertools

#############################################
###
###  Check if a file is not empty
###
#############################################

def fileIsNotEmpty(fpath):
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

#############################################
###
###  Extracting allele id information
###
#############################################

def readIdFile(idFile):
    
    id_list = []
    
    with open(idFile, "r") as input_id:
        # read file containing locus allele IDs.
        for input_id_line in input_id:
            each_input_id = input_id_line.strip("\n")
            id_list.append(each_input_id)
    
    return id_list

#############################################
###
###  Extracting allele seq information
###
#############################################

def readSeqFile(SeqFile):
    
    sequence_list = []
    
    with open(SeqFile, "r") as input_sequence:
        for input_sequence_line in input_sequence:
            each_input_sequence = input_sequence_line.strip("\n")
            sequence_list.append(each_input_sequence)
    
    return sequence_list

#############################################
###
###  Extracting exon bondary information
###
#############################################

def readExonBoundary(exonBoundaryFile):
    
    exonBoundaryDict = {}
    
    with open(exonBoundaryFile, "r") as input_exons_boundary:
        
        for input_exons_line in input_exons_boundary:
            
            each_template         = re.split("\s+", input_exons_line.strip("\n")) 
            exons_boundary_info   = each_template[1]
            exons_positions       = exons_boundary_info.split(",")
            exon2                 = [exons_positions[1], int(exons_positions[2]) - int(exons_positions[1]), exons_positions[2]]
            exon3                 = [exons_positions[2], int(exons_positions[3]) - int(exons_positions[2]), exons_positions[3]]
            exon4                 = [exons_positions[3], int(exons_positions[4]) - int(exons_positions[3]), exons_positions[4]]
            
            exonBoundaryDict[each_template[0]] = [exon2, exon3, exon4]
    
    return exonBoundaryDict

#############################################
###
###  Check if sam file exists
###
#############################################

def checkSAMFile(scratch_folder, output_id, equipment_id):
    
    SAM_file = scratch_folder + "/map_single/" + output_id + ".report.bowtie2.sam"
    SAM_header = []
    
    if not os.path.exists(SAM_file):
        return [False, SAM_header]
    else:
        with open(SAM_file, "r") as f:
            for line in f:
                if equipment_id in line:
                    break
                SAM_header.append(line)     
        return [True, SAM_header]
        
#############################################
###
###  Parsing SAM file
###
#############################################

def parseSAM(sequence):
    
    alignment_detail_list = re.split("\s+", sequence.strip("\n"))
    alignment_start       = alignment_detail_list[3]
    alignment_CIGAR_list  = alignment_detail_list[5].split("M")
    alignment_length      = 0
    
    frontSoftClipped      = 1 #Assuming it has front soft clipped
    endSoftClipped        = 1 #Assuming it has end soft clipped
                                     
    if ("S" in alignment_CIGAR_list[0]) and ("S" in alignment_CIGAR_list[1]):
        
        frontSoftClipped_TRUE = alignment_CIGAR_list[0].split("S")
        endSoftClipped_TRUE   = alignment_CIGAR_list[1].split("S")
        
        frontSoftClipped = frontSoftClipped_TRUE[0]   
        alignment_length = frontSoftClipped_TRUE[1]
        endSoftClipped   = endSoftClipped_TRUE[0]
    
    if ("S" in alignment_CIGAR_list[0]) and ("S" not in alignment_CIGAR_list[1]):
        
        frontSoftClipped_TRUE = alignment_CIGAR_list[0].split("S")
        frontSoftClipped = frontSoftClipped_TRUE[0]   
        alignment_length = frontSoftClipped_TRUE[1]
        endSoftClipped   = 0
    
    if ("S" not in alignment_CIGAR_list[0]) and ("S" in alignment_CIGAR_list[1]):
        
        endSoftClipped_TRUE   = alignment_CIGAR_list[1].split("S")
        
        frontSoftClipped = 0  
        alignment_length = alignment_CIGAR_list[0]
        endSoftClipped   = endSoftClipped_TRUE[0]
    
    if ("S" not in alignment_CIGAR_list[0]) and ("S" not in alignment_CIGAR_list[1]):
        
        frontSoftClipped = 0   
        alignment_length = alignment_CIGAR_list[0]
        endSoftClipped   = 0
                    
    alignment_end = int(alignment_start) + int(alignment_length) - 1
                    
    return [int(alignment_start), int(alignment_end), alignment_detail_list[0], alignment_detail_list[9], int(frontSoftClipped), int(endSoftClipped)]

#############################################
###
###  Calculate Alignment
###
#############################################

def calculateAlignmentLength(alignmentStart, alignmentEnd, exonStart, exonEnd):

    crossingBoundaryFlag = False
    alignedOnExonFlag    = False
    
    alignmentStart = int(alignmentStart)
    alignmentEnd   = int(alignmentEnd)
    exonStart      = int(exonStart)
    exonEnd        = int(exonEnd)
    
    if (alignmentEnd < exonStart or alignmentStart >= exonEnd):
        return [crossingBoundaryFlag, alignedOnExonFlag, int(alignmentStart), int(alignmentEnd) - int(alignmentStart) + 1, int(alignmentEnd), "Case1", 0]
    
    
    elif (alignmentStart >= exonStart and alignmentEnd < exonEnd):
        alignedOnExonFlag = True
        return [crossingBoundaryFlag, alignedOnExonFlag, int(alignmentStart), int(alignmentEnd) - int(alignmentStart) + 1, int(alignmentEnd), "Case2", 0]
    
    
    elif (alignmentStart < exonStart and alignmentEnd >= exonStart and alignmentEnd < exonEnd):
        crossingBoundaryFlag = True
        if (alignmentEnd - exonStart) > (exonStart - alignmentStart):
            alignedOnExonFlag = True
            return [crossingBoundaryFlag, alignedOnExonFlag, int(exonStart), int(alignmentEnd) - int(exonStart) + 1, int(alignmentEnd), "Case3.1", exonStart - alignmentStart + 1]
        else:
            return [crossingBoundaryFlag, alignedOnExonFlag, int(alignmentStart), int(exonStart) - int(alignmentStart) + 1, int(exonStart) - 1, "Case3.2", int(alignmentEnd) - int(exonStart) + 1]
    
    
    elif (alignmentStart >= exonStart and alignmentStart < exonEnd and alignmentEnd >= exonEnd):
        crossingBoundaryFlag = True
        if (exonEnd - alignmentStart) > (alignmentEnd - exonEnd):
            alignedOnExonFlag = True
            return [crossingBoundaryFlag, alignedOnExonFlag, int(alignmentStart), int(exonEnd) - int(alignmentStart), int(exonEnd) - 1, "Case4.1", alignmentEnd - exonEnd + 1]
        else:
            return [crossingBoundaryFlag, alignedOnExonFlag, int(exonEnd), int(alignmentEnd) - int(exonEnd) + 1, int(alignmentEnd), "Case4.2", int(exonEnd) - int(alignmentStart)]
    
    
    elif (alignmentStart < exonStart and exonEnd <= alignmentEnd):
        crossingBoundaryFlag = True
        alignedOnExonFlag = True
        return [crossingBoundaryFlag, alignedOnExonFlag, int(exonStart), int(exonEnd) - int(exonStart), int(exonEnd) - 1, "Case5", 1]
    
    
    else:
        return [crossingBoundaryFlag, alignedOnExonFlag, int(alignmentStart), int(alignmentEnd) - int(alignmentStart) + 1, int(alignmentEnd), "Case6", 0]
        
#############################################
###
###  Calculate Depth of Coverage
###
#############################################

#depthOfCoverage_dict = calculateDepthOfCoverage(exons_boundary_dict, effective_mapped_reads_dict, bothPairReadsMapped_dict, project_folder)

def calculateDepthOfCoverage(exons_boundary_dict, effective_mapped_reads_dict, ID_dictionary, scratch_folder, SAM_header_dict):
    
    DOC_dictionary             = {}
    filtered_mapped_reads_dict = {}
    
    for id in ID_dictionary:
        
        output = open(scratch_folder + "/SAM_FILES/" + id + "_paired.Mapping.sam", "w")
        output.write("".join(SAM_header_dict[id]))
        
        exon2start  = int(exons_boundary_dict[id][0][0])
        exon3start  = int(exons_boundary_dict[id][1][0])
        exon4start  = int(exons_boundary_dict[id][2][0])
        
        exon2length = int(exons_boundary_dict[id][0][1])
        exon3length = int(exons_boundary_dict[id][1][1])
        exon4length = int(exons_boundary_dict[id][2][1])
    
        exon2end    = int(exons_boundary_dict[id][0][2])
        exon3end    = int(exons_boundary_dict[id][1][2])
        exon4end    = int(exons_boundary_dict[id][2][2])
        
        exon2DOC_dictionary = {}
        exon3DOC_dictionary = {}
        exon4DOC_dictionary = {}
        
        exon2_filter_dict   = {}
        exon3_filter_dict   = {}
        exon4_filter_dict   = {}
        
        if (exon2length > 0):
            for key in range(exon2start, exon2end):
                exon2DOC_dictionary[key] = 0
        else:
            exon2DOC_dictionary = { exon2start:0 }
        
        if (exon3length > 0):
            for key in range(exon3start, exon3end):
                exon3DOC_dictionary[key] = 0
        else:
            exon3DOC_dictionary = { exon3start:0 }
        
        if (exon4length > 0):
             for key in range(exon4start, exon4end):
                exon4DOC_dictionary[key] = 0
        else:
            exon4DOC_dictionary = { exon4start:0 }

        if len(ID_dictionary[id]["exon2"]):
            
            for read in ID_dictionary[id]["exon2"]:
                
                if (len(ID_dictionary[id]["exon2"][read]) == 2):
                    for base1 in range(ID_dictionary[id]["exon2"][read][0][0], ID_dictionary[id]["exon2"][read][0][2] + 1):
                        exon2DOC_dictionary[base1] += 1
                    for base2 in range(ID_dictionary[id]["exon2"][read][1][0], ID_dictionary[id]["exon2"][read][1][2] + 1):
                        exon2DOC_dictionary[base2] += 1
                    output.write(ID_dictionary[id]["exon2"][read][0][3])
                    output.write(ID_dictionary[id]["exon2"][read][1][3])
                    
                    exon2_filter_dict[read + ID_dictionary[id]["exon2"][read][0][4]] = [ID_dictionary[id]["exon2"][read][0][0], ID_dictionary[id]["exon2"][read][0][1], ID_dictionary[id]["exon2"][read][0][2]]
                    exon2_filter_dict[read + ID_dictionary[id]["exon2"][read][1][4]] = [ID_dictionary[id]["exon2"][read][1][0], ID_dictionary[id]["exon2"][read][1][1], ID_dictionary[id]["exon2"][read][1][2]]
                    
        else:
            exon2_filter_dict = ID_dictionary[id]["exon2"]
                #else:
                    #output.write(id + " Exon2: " + read + "\n")
                    
        
        if len(ID_dictionary[id]["exon3"]):
            
            for read in ID_dictionary[id]["exon3"]:
                
                if (len(ID_dictionary[id]["exon3"][read]) == 2):
                    for base1 in range(ID_dictionary[id]["exon3"][read][0][0], ID_dictionary[id]["exon3"][read][0][2] + 1):
                        exon3DOC_dictionary[base1] += 1
                    for base2 in range(ID_dictionary[id]["exon3"][read][1][0], ID_dictionary[id]["exon3"][read][1][2] + 1):
                        exon3DOC_dictionary[base2] += 1
                    output.write(ID_dictionary[id]["exon3"][read][0][3])
                    output.write(ID_dictionary[id]["exon3"][read][1][3])
                    
                    exon3_filter_dict[read + ID_dictionary[id]["exon3"][read][0][4]] = [ID_dictionary[id]["exon3"][read][0][0], ID_dictionary[id]["exon3"][read][0][1], ID_dictionary[id]["exon3"][read][0][2]]
                    exon3_filter_dict[read + ID_dictionary[id]["exon3"][read][1][4]] = [ID_dictionary[id]["exon3"][read][1][0], ID_dictionary[id]["exon3"][read][1][1], ID_dictionary[id]["exon3"][read][1][2]]
        
        else:
            exon3_filter_dict = ID_dictionary[id]["exon3"]
                #else:
                    #output.write(id + " Exon3: " + read + "\n")
               
        
        if len(ID_dictionary[id]["exon4"]):
            
            for read in ID_dictionary[id]["exon4"]:
                
                if (len(ID_dictionary[id]["exon4"][read]) == 2):
                    for base1 in range(ID_dictionary[id]["exon4"][read][0][0], ID_dictionary[id]["exon4"][read][0][2] + 1):
                        exon4DOC_dictionary[base1] += 1
                    for base2 in range(ID_dictionary[id]["exon4"][read][1][0], ID_dictionary[id]["exon4"][read][1][2] + 1):
                        exon4DOC_dictionary[base2] += 1
                    output.write(ID_dictionary[id]["exon4"][read][0][3])
                    output.write(ID_dictionary[id]["exon4"][read][1][3])
                    
                    exon4_filter_dict[read + ID_dictionary[id]["exon4"][read][0][4]] = [ID_dictionary[id]["exon4"][read][0][0], ID_dictionary[id]["exon4"][read][0][1], ID_dictionary[id]["exon4"][read][0][2]]
                    exon4_filter_dict[read + ID_dictionary[id]["exon4"][read][1][4]] = [ID_dictionary[id]["exon4"][read][1][0], ID_dictionary[id]["exon4"][read][1][1], ID_dictionary[id]["exon4"][read][1][2]]
        
        else:
            exon4_filter_dict = ID_dictionary[id]["exon4"]
          
        DOC_dictionary[id] = {
                              "exon2" : exon2DOC_dictionary,
                              "exon3" : exon3DOC_dictionary,
                              "exon4" : exon4DOC_dictionary
                              }
        
        filtered_mapped_reads_dict[id] = {
                                          "exon2" : exon2_filter_dict,
                                          "exon3" : exon3_filter_dict,
                                          "exon4" : exon4_filter_dict
                                          }
    
        output.close() 
        
    return [DOC_dictionary, filtered_mapped_reads_dict]

#############################################
###
###  Print out depth of coverage
###
#############################################

def printOutDepthOfCoverage(depthOfCoverage_dict, project_folder):
    
    output = open(project_folder + "/detail/detail_depthOfCoverage.txt", "w")
    
    for id in depthOfCoverage_dict:
        
        exon2Positions = []
        exon2Coverages = []
        exon3Positions = []
        exon3Coverages = []
        exon4Positions = []
        exon4Coverages = []
        
        exon2Positions = [id + "_Exon2"] + [position for position in sorted(depthOfCoverage_dict[id]["exon2"].iterkeys())]
        exon3Positions = [id + "_Exon3"] + [position for position in sorted(depthOfCoverage_dict[id]["exon3"].iterkeys())]
        exon4Positions = [id + "_Exon4"] + [position for position in sorted(depthOfCoverage_dict[id]["exon4"].iterkeys())]
        
        exon2Coverages = [id + "_Exon2"] + [depthOfCoverage_dict[id]["exon2"][position] for position in sorted(depthOfCoverage_dict[id]["exon2"].iterkeys())]
        exon3Coverages = [id + "_Exon3"] + [depthOfCoverage_dict[id]["exon3"][position] for position in sorted(depthOfCoverage_dict[id]["exon3"].iterkeys())]
        exon4Coverages = [id + "_Exon4"] + [depthOfCoverage_dict[id]["exon4"][position] for position in sorted(depthOfCoverage_dict[id]["exon4"].iterkeys())]
        
        output.write("\t".join([str(i) for i in exon2Positions]) + "\n")
        output.write("\t".join([str(i) for i in exon2Coverages]) + "\n")
        
        output.write("\t".join([str(i) for i in exon3Positions]) + "\n")
        output.write("\t".join([str(i) for i in exon3Coverages]) + "\n")
        
        output.write("\t".join([str(i) for i in exon4Positions]) + "\n")
        output.write("\t".join([str(i) for i in exon4Coverages]) + "\n")
    
    output.close()
    
    if fileIsNotEmpty(project_folder + "/detail/detail_depthOfCoverage.txt"):
        return True
    else:
        return False
    
#############################################
###
###  Determine Dynamic Boundary
###
#############################################

def dynamicBoundary(Coverage):
    
    if len(Coverage) <= 1:
        #print "Exon missing."
        return [0,0]
    
    if max(Coverage) == 0:
        #print "No reads mapped on this exon."
        return [0,0]
    
    left_boundary  = 0
    right_coundary = 0
    
    #Coverage       = Coverage[10:-10]
    reverseList    = Coverage[::-1]
    
    left_boundary  = Coverage.index( filter(lambda x: x!=0, Coverage)[0] )
    right_coundary = reverseList.index( filter(lambda x: x!=0, reverseList)[0] )
        
        #print "Left boundary: ", left_boundary
        #print "Right boundary: ", right_coundary
        
    if left_boundary > 0 and right_coundary > 0:
        return [min(Coverage[left_boundary:-right_coundary]), 0]
    elif left_boundary == 0 and right_coundary > 0:
        return [min(Coverage[:-right_coundary]), 0]
    elif left_boundary > 0 and right_coundary == 0:
        return [min(Coverage[left_boundary:]), 0]
    elif left_boundary == 0 and right_coundary == 0:
        return [min(Coverage), 1]
    else:
            #print "Error: can not determin the boundary!"
        return [0,0]

#############################################
###
###  Compute minimum coverage
###
#############################################

def minimumCoverage(depthOfCoverage_dict):
    
    minimumCoverage_dictionary = {}
    
    for id in depthOfCoverage_dict:
        
        exon2Coverage = [depthOfCoverage_dict[id]["exon2"][base] for base in sorted(depthOfCoverage_dict[id]["exon2"].iterkeys())]
        exon3Coverage = [depthOfCoverage_dict[id]["exon3"][base] for base in sorted(depthOfCoverage_dict[id]["exon3"].iterkeys())]
        exon4Coverage = [depthOfCoverage_dict[id]["exon4"][base] for base in sorted(depthOfCoverage_dict[id]["exon4"].iterkeys())]
        
        #print id
        #print exon2Coverage
        #print exon3Coverage
        #print exon4Coverage
        
        exon2MinimumCoverage = dynamicBoundary(exon2Coverage)
        exon3MinimumCoverage = dynamicBoundary(exon3Coverage)
        exon4MinimumCoverage = dynamicBoundary(exon4Coverage)
        
        #print id
        #print exon2MinimumCoverage
        #print exon3MinimumCoverage
        #print exon4MinimumCoverage
        
        minimumCoverage_dictionary[id] = {
                                          "exon2" : exon2MinimumCoverage,
                                          "exon3" : exon3MinimumCoverage,
                                          "exon4" : exon4MinimumCoverage
                                          }
    return minimumCoverage_dictionary

#############################################
###
###  Print out whole exon coverage
###
#############################################

def printOutWholeExonCoverage(exonWholeCoverage_dict, project_folder):
    
    output = open(project_folder + "/detail/detail_WholeExonCoverage.txt", "w")
    output.write("Allele\tExon2\tExon3\tExon4\n")
    
    for id in exonWholeCoverage_dict:
        
        output.write(id + "\t")
        output.write(str(exonWholeCoverage_dict[id][0]) + "\t")
        output.write(str(exonWholeCoverage_dict[id][1]) + "\t")
        output.write(str(exonWholeCoverage_dict[id][2]) + "\n")
    
    output.close()
    
    if fileIsNotEmpty(project_folder + "/detail/detail_WholeExonCoverage.txt"):
        return True
    else:
        return False
 

#############################################
###
###  Print out minimum coverage
###
#############################################

def printOutMinimumCoverage(minimumCoverage_dict, project_folder):
    
    output = open(project_folder + "/detail/detail_minimumCoverage.txt", "w")
    output.write("Allele\tExon2\tExon2_Detail\tExon3\tExon3_Detail\tExon4\tExon4_Detail\n")
    
    for id in minimumCoverage_dict:
        
        output.write(id + "\t")
        output.write(str(minimumCoverage_dict[id]["exon2"][0]) + "\t")
        output.write(str(minimumCoverage_dict[id]["exon2"][1]) + "\t")
        output.write(str(minimumCoverage_dict[id]["exon3"][0]) + "\t")
        output.write(str(minimumCoverage_dict[id]["exon3"][1]) + "\t")
        output.write(str(minimumCoverage_dict[id]["exon4"][0]) + "\t")
        output.write(str(minimumCoverage_dict[id]["exon4"][1]) + "\n")
    
    output.close()
    
    if fileIsNotEmpty(project_folder + "/detail/detail_minimumCoverage.txt"):
        return True
    else:
        return False
        
#############################################
###
###  Generate alignment info
###
#############################################

def generateExonAlignmentInfo(scratch_folder, output_id, exons_boundary_dict, equipment_id, alignmentLengthCutoff):
    
    AlignmentDict            = {}
    crossingExonBoundary     = {}
    exon2Alignment           = {}
    exon3Alignment           = {}   
    exon4Alignment           = {}
    
    pairedReads_mapping_dict = {}
    exon2_pairedReads        = {}
    exon3_pairedReads        = {}
    exon4_pairedReads        = {}
    
    exon2WholeCoverage   = 0
    exon3WholeCoverage   = 0
    exon4WholeCoverage   = 0
    
    exon2start  = int(exons_boundary_dict[output_id][0][0])
    exon3start  = int(exons_boundary_dict[output_id][1][0])
    exon4start  = int(exons_boundary_dict[output_id][2][0])
        
    exon2length = int(exons_boundary_dict[output_id][0][1])
    exon3length = int(exons_boundary_dict[output_id][1][1])
    exon4length = int(exons_boundary_dict[output_id][2][1])
    
    exon2end    = int(exons_boundary_dict[output_id][0][2])
    exon3end    = int(exons_boundary_dict[output_id][1][2])
    exon4end    = int(exons_boundary_dict[output_id][2][2])
    
    SAM_file      = scratch_folder + "/map_single/" + output_id + ".report.bowtie2.sam"
    #test_SAM_file = open(scratch_folder + "/" + output_id + ".report.bowtie2.sam", "w")
    
    with open(SAM_file, "r") as alignment_results:
        
        for each_aligned_sequence in alignment_results:
            
            alignment_start   = 0
            alignment_end     = 0
            shortReadID       = ""
            shortReadSeq      = "" 
            frontSoftClipped  = 1
            endSoftClipped    = 1
            
            matched_align_info = re.search(equipment_id, each_aligned_sequence)
            
            if matched_align_info is not None:     
               [alignment_start, alignment_end, shortReadID, shortReadSeq, frontSoftClipped, endSoftClipped] = parseSAM(each_aligned_sequence)
            else:
                continue
            
            if (exon2length > 0):
                
                exon2_List = calculateAlignmentLength(alignment_start, alignment_end, exon2start, exon2end)
                
                if (exon2_List[1] and exon2_List[3] >= alignmentLengthCutoff):
                    
                    if (exon2start < exon2_List[2]) and (exon2_List[4] < exon2end - 1) and (not frontSoftClipped) and (not endSoftClipped):
                        exon2Alignment[shortReadID + "|" + shortReadSeq] = [exon2_List[2], exon2_List[3], exon2_List[4]]
                        #test_SAM_file.write("Exon2 : *Case1 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon2_pairedReads:
                            exon2_pairedReads[shortReadID] = []
                            exon2_pairedReads[shortReadID].append([exon2_List[2], exon2_List[3], exon2_List[4], each_aligned_sequence, shortReadSeq])
                        else:
                            exon2_pairedReads[shortReadID].append([exon2_List[2], exon2_List[3], exon2_List[4], each_aligned_sequence, shortReadSeq])
                        
                    
                    if (exon2start < exon2_List[2]) and (exon2_List[4] == exon2end - 1) and (not frontSoftClipped):
                        exon2Alignment[shortReadID + "|" + shortReadSeq] = [exon2_List[2], exon2_List[3], exon2_List[4]]
                        #test_SAM_file.write("Exon2 : *Case2 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon2_pairedReads:
                            exon2_pairedReads[shortReadID] = []
                            exon2_pairedReads[shortReadID].append([exon2_List[2], exon2_List[3], exon2_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon2_pairedReads[shortReadID].append([exon2_List[2], exon2_List[3], exon2_List[4],each_aligned_sequence, shortReadSeq])
                            
                    
                    if (exon2start == exon2_List[2]) and (exon2_List[4] < exon2end - 1) and (not endSoftClipped):
                        exon2Alignment[shortReadID + "|" + shortReadSeq] = [exon2_List[2], exon2_List[3], exon2_List[4]]
                        #test_SAM_file.write("Exon2 : *Case3 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon2_pairedReads:
                            exon2_pairedReads[shortReadID] = []
                            exon2_pairedReads[shortReadID].append([exon2_List[2], exon2_List[3], exon2_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon2_pairedReads[shortReadID].append([exon2_List[2], exon2_List[3], exon2_List[4],each_aligned_sequence, shortReadSeq])
                    
                    if (exon2start == exon2_List[2]) and (exon2_List[4] == exon2end - 1):
                        exon2Alignment[shortReadID + "|" + shortReadSeq] = [exon2_List[2], exon2_List[3], exon2_List[4]]
                        exon2WholeCoverage += 1
                        #test_SAM_file.write("Exon2 : *Case4 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon2_pairedReads:
                            exon2_pairedReads[shortReadID] = []
                            exon2_pairedReads[shortReadID].append([exon2_List[2], exon2_List[3], exon2_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon2_pairedReads[shortReadID].append([exon2_List[2], exon2_List[3], exon2_List[4],each_aligned_sequence, shortReadSeq])
                        
                if (exon2_List[0]):
                    crossingExonBoundary["exon2 " + each_aligned_sequence] = exon2_List + ["exon2"]
                    
            
            if (exon3length > 0):
                
                exon3_List = calculateAlignmentLength(alignment_start, alignment_end, exon3start, exon3end)
                
                if (exon3_List[1] and exon3_List[3] >= alignmentLengthCutoff):
                    
                    if (exon3start < exon3_List[2]) and (exon3_List[4] < exon3end - 1) and (not frontSoftClipped) and (not endSoftClipped):
                        exon3Alignment[shortReadID + "|" + shortReadSeq] = [exon3_List[2], exon3_List[3], exon3_List[4]]
                        #test_SAM_file.write("Exon3 : *Case1 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon3_pairedReads:
                            exon3_pairedReads[shortReadID] = []
                            exon3_pairedReads[shortReadID].append([exon3_List[2], exon3_List[3], exon3_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon3_pairedReads[shortReadID].append([exon3_List[2], exon3_List[3], exon3_List[4],each_aligned_sequence, shortReadSeq])
                    
                    if (exon3start < exon3_List[2]) and (exon3_List[4] == exon3end - 1) and (not frontSoftClipped):
                        exon3Alignment[shortReadID + "|" + shortReadSeq] = [exon3_List[2], exon3_List[3], exon3_List[4]]
                        #test_SAM_file.write("Exon3 : *Case2 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon3_pairedReads:
                            exon3_pairedReads[shortReadID] = []
                            exon3_pairedReads[shortReadID].append([exon3_List[2], exon3_List[3], exon3_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon3_pairedReads[shortReadID].append([exon3_List[2], exon3_List[3], exon3_List[4],each_aligned_sequence, shortReadSeq])
                    
                    if (exon3start == exon3_List[2]) and (exon3_List[4] < exon3end - 1) and (not endSoftClipped):
                        exon3Alignment[shortReadID + "|" + shortReadSeq] = [exon3_List[2], exon3_List[3], exon3_List[4]]
                        #test_SAM_file.write("Exon3 : *Case3 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon3_pairedReads:
                            exon3_pairedReads[shortReadID] = []
                            exon3_pairedReads[shortReadID].append([exon3_List[2], exon3_List[3], exon3_List[4], each_aligned_sequence, shortReadSeq])
                        else:
                            exon3_pairedReads[shortReadID].append([exon3_List[2], exon3_List[3], exon3_List[4], each_aligned_sequence, shortReadSeq])
                    
                    
                    if (exon3start == exon3_List[2]) and (exon3_List[4] == exon3end - 1):
                        exon3Alignment[shortReadID + "|" + shortReadSeq] = [exon3_List[2], exon3_List[3], exon3_List[4]]
                        exon3WholeCoverage += 1
                        #test_SAM_file.write("Exon3 : *Case4 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon3_pairedReads:
                            exon3_pairedReads[shortReadID] = []
                            exon3_pairedReads[shortReadID].append([exon3_List[2], exon3_List[3], exon3_List[4] ,each_aligned_sequence, shortReadSeq])
                        else:
                            exon3_pairedReads[shortReadID].append([exon3_List[2], exon3_List[3], exon3_List[4] ,each_aligned_sequence, shortReadSeq])
                    
                
                if (exon3_List[0]):
                    crossingExonBoundary["exon3 " + each_aligned_sequence] = exon3_List + ["exon3"]
            
            if (exon4length > 0):
                
                exon4_List = calculateAlignmentLength(alignment_start, alignment_end, exon4start, exon4end)
                
                if (exon4_List[1] and exon4_List[3] >= alignmentLengthCutoff):
                    
                    if (exon4start < exon4_List[2]) and (exon4_List[4] < exon4end - 1) and (not frontSoftClipped) and (not endSoftClipped):
                        exon4Alignment[shortReadID + "|" + shortReadSeq] = [exon4_List[2], exon4_List[3], exon4_List[4]]
                        #test_SAM_file.write("Exon4 : *Case1 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon4_pairedReads:
                            exon4_pairedReads[shortReadID] = []
                            exon4_pairedReads[shortReadID].append([exon4_List[2], exon4_List[3], exon4_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon4_pairedReads[shortReadID].append([exon4_List[2], exon4_List[3], exon4_List[4],each_aligned_sequence, shortReadSeq])
                    
                    if (exon4start < exon4_List[2]) and (exon4_List[4] == exon4end - 1) and (not frontSoftClipped):
                        exon4Alignment[shortReadID + "|" + shortReadSeq] = [exon4_List[2], exon4_List[3], exon4_List[4]]
                        #test_SAM_file.write("Exon4 : *Case2 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon4_pairedReads:
                            exon4_pairedReads[shortReadID] = []
                            exon4_pairedReads[shortReadID].append([exon4_List[2], exon4_List[3], exon4_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon4_pairedReads[shortReadID].append([exon4_List[2], exon4_List[3], exon4_List[4],each_aligned_sequence, shortReadSeq])
                    
                    
                    if (exon4start == exon4_List[2]) and (exon4_List[4] < exon4end - 1) and (not endSoftClipped):
                        exon4Alignment[shortReadID + "|" + shortReadSeq] = [exon4_List[2], exon4_List[3], exon4_List[4]]
                        #test_SAM_file.write("Exon4 : *Case3 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon4_pairedReads:
                            exon4_pairedReads[shortReadID] = []
                            exon4_pairedReads[shortReadID].append([exon4_List[2], exon4_List[3], exon4_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon4_pairedReads[shortReadID].append([exon4_List[2], exon4_List[3], exon4_List[4],each_aligned_sequence, shortReadSeq])
                    
                    
                    if (exon4start == exon4_List[2]) and (exon4_List[4] == exon4end - 1):
                        exon4Alignment[shortReadID + "|" + shortReadSeq] = [exon4_List[2], exon4_List[3], exon4_List[4]]
                        exon4WholeCoverage += 1
                        #test_SAM_file.write("Exon4 : *Case4 \n" +  each_aligned_sequence)
                        
                        if shortReadID not in exon4_pairedReads:
                            exon4_pairedReads[shortReadID] = []
                            exon4_pairedReads[shortReadID].append([exon4_List[2], exon4_List[3], exon4_List[4],each_aligned_sequence, shortReadSeq])
                        else:
                            exon4_pairedReads[shortReadID].append([exon4_List[2], exon4_List[3], exon4_List[4],each_aligned_sequence, shortReadSeq])
                    
                
                if (exon4_List[0]):
                    crossingExonBoundary["exon4 " + each_aligned_sequence] = exon4_List + ["exon4"]
    
    AlignmentDict = {
                     "exon2" : exon2Alignment,
                     "exon3" : exon3Alignment,
                     "exon4" : exon4Alignment,
                     "crossingExonBoundary" : crossingExonBoundary
                     }
    
    pairedReads_mapping_dict = {
                                "exon2" : exon2_pairedReads,
                                "exon3" : exon3_pairedReads,
                                "exon4" : exon4_pairedReads
                                }
    
    #test_SAM_file.close()
    return [AlignmentDict, [exon2WholeCoverage, exon3WholeCoverage, exon4WholeCoverage], pairedReads_mapping_dict]

################################################################
#
#  Print out reads crossing the boundary
#
################################################################

def printOutcrossingExonBoundary(effective_mapped_reads_dict, project_folder):
    
    output = open(project_folder + "/detail/detail_crossingExonBoundary.txt", "w")
    
    for id in effective_mapped_reads_dict:
        
        for read in effective_mapped_reads_dict[id]["crossingExonBoundary"]:
            output.write(id + ":  ")
            output.write(read + "\n")
            output.write(str(effective_mapped_reads_dict[id]["crossingExonBoundary"][read]) + "\n")
            
    output.close()

################################################################
#
#  Print out select cDNA reference
#
################################################################
    
def printOutAlleleList(cDNA_list, project_folder):
    
    output = open(project_folder + "/remain_id.txt", "w")
    
    for cDNA in cDNA_list:
        
        output.write(str(cDNA[0]) + "\t" + str(cDNA[1]) + "\n")
    
    output.close()
    
################################################################
#
#  Print out select cDNA reference
#
################################################################

def generateBamFiles(cDNA_list, project_folder, scratch_folder):
    
    samtools         = "/share/apps/samtools/samtools"
    BAM_File_counter = 0
    
    for id in cDNA_list:
        
        input = scratch_folder + "/SAM_FILES/" + id + "_paired.Mapping.sam"
        output = project_folder + "/SAM_FILES/" + id + "_paired.Mapping"

        #samtools view -bS DPB1_13\:01.report.bowtie2.sam | samtools sort - DPB1_13_01
        #samtools index DPB1_13_01.bam DPB1_13_01.bai
        
        command1 = samtools + " view -bS " + input + " | " + samtools + " sort - " + output
        command2 = samtools + " index " + output + ".bam " + output + ".bai"
        
        os.system(command1)
        os.system(command2)
        
        if fileIsNotEmpty(output + ".bam"):
            BAM_File_counter += 1
    
    if (len(cDNA_list) == BAM_File_counter):
        return True
    else:
        return False
        
################################################################
#
#  Select cDNA reference
#
################################################################

def selectSingleAllele(minimumCoverage, exonWholeCoverage_dict, project_folder, primer_start, primer_end):
    
    print "Start to select single cDNA.\n"
    print "Primer Start: ", primer_start
    print "Primer end  : ", primer_end
    
    cDNA_list             = []
    selectedNumberCutoff  = 0
    selectedID_list       = []
    fourDigitsAllele_List = []
    
    if primer_start == 2 and primer_end == 4:
    
        for id in exonWholeCoverage_dict:
            
            exonWholeCoverage_exon2 = int(exonWholeCoverage_dict[id][0])
            exonWholeCoverage_exon3 = int(exonWholeCoverage_dict[id][1])
            exonWholeCoverage_exon4 = int(exonWholeCoverage_dict[id][2])
            exonWholeCoverage = exonWholeCoverage_exon2 * exonWholeCoverage_exon3 * exonWholeCoverage_exon4
            
            if (exonWholeCoverage_exon2 > 1 and exonWholeCoverage_exon3 > 1 and exonWholeCoverage_exon4 > 1):
                cDNA_list.append([id, exonWholeCoverage])
                selectedID_list.append(id)
                [loci, alleleName] = id.split("_")
                alleleName_list = alleleName.split(":")
                allele4digits = alleleName_list[0] + ":" + alleleName_list[1]
                print "4 digits:", allele4digits
                fourDigitsAllele_List.append(allele4digits)
                print id + "\t" + str(exonWholeCoverage) + "\t" + str(exonWholeCoverage_dict[id][0]) + "\t" + str(exonWholeCoverage_dict[id][1]) + "\t" + str(exonWholeCoverage_dict[id][2]) + "\n"
            else:
                continue
            
        print "We have %d alleles that have single read covering the whole exon on all three exons."%len(cDNA_list)
        fourDigitsAllele_Set = set(fourDigitsAllele_List)
        print "We have %d alleles in 4 digits allele set."%len(fourDigitsAllele_Set)
        
        if len(fourDigitsAllele_Set) < 3:
        
            for id in minimumCoverage:
        
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0], minimumCoverage[id]["exon4"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
                coveredWholeExon  = int(minimumCoverage[id]["exon2"][1]) + int(minimumCoverage[id]["exon3"][1]) + int(minimumCoverage[id]["exon4"][1])
        
                if count_minCoverage == 3 and coveredWholeExon == 3 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                    print id + str(minCoverageList) + "\n"
                else:
                    continue
    
        print "We have %d alleles that all three exons were totally covered by reads."%len(cDNA_list)
        
    
        if (len(cDNA_list) == 0):
        
            for id in minimumCoverage:
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0], minimumCoverage[id]["exon4"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
             
                if (count_minCoverage == 3) and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that have mappability on all three exons."%len(cDNA_list)
    
        if (len(cDNA_list) == 0):
        
            for id in minimumCoverage:
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0], minimumCoverage[id]["exon4"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
                coveredWholeExon  = int(minimumCoverage[id]["exon2"][1]) + int(minimumCoverage[id]["exon3"][1]) + int(minimumCoverage[id]["exon4"][1])
        
                if count_minCoverage == 2 and coveredWholeExon == 2 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that two exons were totally covered by reads."%len(cDNA_list)
    
        if (len(cDNA_list) == 0):
        
            for id in minimumCoverage:
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0], minimumCoverage[id]["exon4"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
         
                if count_minCoverage == 2 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that have mappability on two exons."%len(cDNA_list)
    
        if (len(cDNA_list) == 0):
            
            for id in minimumCoverage:
            
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0], minimumCoverage[id]["exon4"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
                coveredWholeExon  = int(minimumCoverage[id]["exon2"][1]) + int(minimumCoverage[id]["exon3"][1]) + int(minimumCoverage[id]["exon4"][1])
        
                if count_minCoverage == 1 and coveredWholeExon == 1 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that one exon were totally covered by reads."%len(cDNA_list)
    
        if (len(cDNA_list) == 0):
            
            for id in minimumCoverage:
            
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0], minimumCoverage[id]["exon4"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
            
                if count_minCoverage == 1 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that have mappability on one exon."%len(cDNA_list)
        
    
    if primer_start == 2 and primer_end == 3:
        
        if (len(cDNA_list) == 0):
        
            for id in minimumCoverage:
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
                coveredWholeExon  = int(minimumCoverage[id]["exon2"][1]) + int(minimumCoverage[id]["exon3"][1])
        
                if count_minCoverage == 2 and coveredWholeExon == 2 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that two exons were totally covered by reads."%len(cDNA_list)
    
        if (len(cDNA_list) == 0):
        
            for id in minimumCoverage:
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
         
                if count_minCoverage == 2 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that have mappability on two exons."%len(cDNA_list)
    
        if (len(cDNA_list) == 0):
            
            for id in minimumCoverage:
            
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
                coveredWholeExon  = int(minimumCoverage[id]["exon2"][1]) + int(minimumCoverage[id]["exon3"][1])
                
                if count_minCoverage == 1 and coveredWholeExon == 1 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that one exon were totally covered by reads."%len(cDNA_list)
    
        if (len(cDNA_list) == 0):
            
            for id in minimumCoverage:
            
                minCoverageList   = [minimumCoverage[id]["exon2"][0], minimumCoverage[id]["exon3"][0]]
                count_minCoverage = sum(x > 0 for x in minCoverageList)
            
                if count_minCoverage == 1 and id not in selectedID_list:
                    selectedID_list.append(id)
                    cDNA_list.append([id, max(minCoverageList)])
                else:
                    continue
    
        print "We have %d alleles that have mappability on one exon."%len(cDNA_list)
    
    
    if (len(cDNA_list) == 0):
        return cDNA_list
    
    if (len(cDNA_list) > 200):
        selectedNumberCutoff = 200
    else:
        selectedNumberCutoff = len(cDNA_list)
    
    cDNA_list.sort(key=operator.itemgetter(1),reverse=True)
    
    selectedAlleles = cDNA_list[:selectedNumberCutoff]
    
    printOutAlleleList(cDNA_list, project_folder)
    
    return [allele[0] for allele in selectedAlleles]

################################################################
#
#  Allele Combination Selection
#
################################################################

def choseAlleleCombination(selectedAllele_list, effective_mapped_reads_dict, primer_start, primer_end):
    
    selectedAlleleComboList = []
    totalMappedReadsList    = []
    
    for allele_combination in list(itertools.combinations(selectedAllele_list, 2)):
        
        allele1 = allele_combination[0]
        allele2 = allele_combination[1]
        
        totalMappedReads     = 0
        totalMappedReadsInfo = []
        exonMappedReadsInfo  = []
        
        totalCoverage        = 0
        totalCoverageInfo    = []
        exonCoverageInfo     = []
        
        totalCoverageInfo.append(allele1 + "_" + allele2)
        totalMappedReadsInfo.append(allele1 + "_" + allele2)
        
        for exon in range(primer_start, primer_end + 1):
            
            allele1_coverage_dict = {}
            allele2_coverage_dict = {}
            
            for read1 in effective_mapped_reads_dict[allele1]["exon" + str(exon)]:
                allele1_coverage_dict[read1] = effective_mapped_reads_dict[allele1]["exon" + str(exon)][read1]
                
            for read2 in effective_mapped_reads_dict[allele2]["exon" + str(exon)]:
                allele2_coverage_dict[read2] = effective_mapped_reads_dict[allele2]["exon" + str(exon)][read2]
        
            allele1_set = set(allele1_coverage_dict.keys())
            allele2_set = set(allele2_coverage_dict.keys())
        
            mapped_allele1_only        = len(allele1_set - allele2_set)
            mapped_allele2_only        = len(allele2_set - allele1_set)
            mapped_allele_intersection = len(allele1_set & allele2_set)
            mapped_allele_union        = mapped_allele1_only + mapped_allele2_only + mapped_allele_intersection
            
            totalMappedReads += mapped_allele_union
            exonMappedReadsInfo.extend([mapped_allele1_only, mapped_allele2_only, mapped_allele_intersection, mapped_allele_union])
            
            allele1_only        = 0
            allele2_only        = 0
            allele_union        = 0
            allele_intersection = 0
            
            for id in (allele1_set - allele2_set):
                allele1_only += allele1_coverage_dict[id][1]
        
            for id in (allele2_set - allele1_set):
                allele2_only += allele2_coverage_dict[id][1]
        
            if len(allele1_set & allele2_set) > 0:
                for id in (allele1_set & allele2_set):
                    if allele1_coverage_dict[id][1] > allele2_coverage_dict[id][1]:
                        allele_intersection += allele1_coverage_dict[id][1]
                    else:
                        allele_intersection += allele2_coverage_dict[id][1]
        
            allele_union = allele1_only + allele2_only + allele_intersection
            totalCoverage += allele_union
            exonCoverageInfo.extend([allele1_only, allele2_only, allele_intersection, allele_union])
            
        totalCoverageInfo.append(totalCoverage)
        totalCoverageInfo.extend(exonCoverageInfo)
        selectedAlleleComboList.append(totalCoverageInfo)
        
        totalMappedReadsInfo.append(totalMappedReads)
        totalMappedReadsInfo.extend(exonMappedReadsInfo)
        totalMappedReadsList.append(totalMappedReadsInfo)
        
    for putative_allele in selectedAllele_list:
        
        totalCoverage     = 0
        totalCoverageInfo = []
        exonCoverageInfo  = []
        totalCoverageInfo.append(putative_allele + "_" + putative_allele)
        
        totalMappedReads     = 0
        totalMappedReadsInfo = []
        exonMappedReadsInfo  = []
        totalMappedReadsInfo.append(putative_allele + "_" + putative_allele)
        
        for exon in range(primer_start, primer_end + 1):
            
            allele1_coverage_dict = {}
            
            for read in effective_mapped_reads_dict[putative_allele]["exon" + str(exon)]:
                allele1_coverage_dict[read] = effective_mapped_reads_dict[putative_allele]["exon" + str(exon)][read]
        
            allele1_set = set(allele1_coverage_dict.keys())
            
            mapped_allele1_only        = 0
            mapped_allele2_only        = 0
            mapped_allele_intersection = len(allele1_set)
            mapped_allele_union        = mapped_allele_intersection * 1.05
            
            totalMappedReads += mapped_allele_union
            exonMappedReadsInfo.extend([mapped_allele1_only, mapped_allele2_only, mapped_allele_intersection, mapped_allele_union])
           
            allele1_only        = 0
            allele2_only        = 0
            allele_union        = 0
            allele_intersection = 0
        
            for id in allele1_set:
                allele_intersection += allele1_coverage_dict[id][1]
        
            allele_union = allele_intersection * 1.05
            totalCoverage += allele_union
            exonCoverageInfo.extend([allele1_only, allele2_only, allele_intersection, allele_union])
        
        totalCoverageInfo.append(totalCoverage)
        totalCoverageInfo.extend(exonCoverageInfo)
        selectedAlleleComboList.append(totalCoverageInfo)
        
        totalMappedReadsInfo.append(totalMappedReads)
        totalMappedReadsInfo.extend(exonMappedReadsInfo)
        totalMappedReadsList.append(totalMappedReadsInfo)
            
    return [selectedAlleleComboList, totalMappedReadsList]

################################################################
#
#  Printout Allele Combination Mapped Reads
#
################################################################

def printOutTotalMappedReads(alleleCombination_list, project_folder):
    
    if (len(alleleCombination_list) == 0):
        return False
    
    alleleCombination_list.sort(key=operator.itemgetter(1),reverse=True)
    
    count       = open(project_folder + "/AlleleCombo_MappedReads.txt", "w")
    countDetail = open(project_folder + "/detail/detail_AlleleCombo_MappedReads.txt", "w")
    
    countDetail.write("AlleleCombos\tTotalMappedReads\t")
    countDetail.write("Exon2Allele1_Only\tExon2Allele2_Only\tExon2Intersection\tExon2Union\t")
    countDetail.write("Exon3Allele1_Only\tExon3Allele2_Only\tExon3Intersection\tExon3Union\t")
    countDetail.write("Exon4Allele1_Only\tExon4Allele2_Only\tExon4Intersection\tExon4Union\n")
    
    for allele in alleleCombination_list:
        count.write(allele[0] + "\t" + str(allele[1]) + "\n")
        countDetail.write("\t".join([str(i) for i in allele]) + "\n")
    
    count.close()
    countDetail.close()
    
    if fileIsNotEmpty(project_folder + "/AlleleCombo_MappedReads.txt") and fileIsNotEmpty(project_folder + "/detail/detail_AlleleCombo_MappedReads.txt"):
        return True
    else:
        return False

################################################################
#
#  Printout Allele Combination Selections
#
################################################################

def printOutAlleleCombo(alleleCombination_list, project_folder):
    
    if (len(alleleCombination_list) == 0):
        return False
    
    alleleCombination_list.sort(key=operator.itemgetter(1),reverse=True)
    
    count       = open(project_folder + "/count_2.txt", "w")
    countDetail = open(project_folder + "/detail/detail_count_2.txt", "w")
    
    countDetail.write("AlleleCombos\tTotalCoverage\t")
    countDetail.write("Exon2Allele1_Only\tExon2Allele2_Only\tExon2Intersection\tExon2Union\t")
    countDetail.write("Exon3Allele1_Only\tExon3Allele2_Only\tExon3Intersection\tExon3Union\t")
    countDetail.write("Exon4Allele1_Only\tExon4Allele2_Only\tExon4Intersection\tExon4Union\n")
    
    for allele in alleleCombination_list:
        count.write(allele[0] + "\t" + str(allele[1]) + "\n")
        countDetail.write("\t".join([str(i) for i in allele]) + "\n")
    
    count.close()
    countDetail.close()
    
    if fileIsNotEmpty(project_folder + "/count_2.txt") and fileIsNotEmpty(project_folder + "/detail/detail_count_2.txt"):
        return True
    else:
        return False
    
################################################################
#
#  Printout Mappability Per Exon
#
################################################################

def printOutMappabilityPerExon(effective_mapped_reads_dict, project_folder):
    
    if (len(effective_mapped_reads_dict) == 0):
        return False
    
    countDetail = open(project_folder + "/detail/mappedReadsPerExon.txt", "w")
    baseDetail  = open(project_folder + "/detail/mappedBasePerExon.txt", "w") 
    
    countDetail.write("Allele\tExon2\tExon3\tExon4\n")
    baseDetail.write("Allele\tExon2\tExon3\tExon4\n")
    
    for id in effective_mapped_reads_dict:
        
        countDetail.write(id + "\t")
        baseDetail.write(id + "\t")
        
        exon2MappedReads = str(len(effective_mapped_reads_dict[id]["exon2"]))
        exon3MappedReads = str(len(effective_mapped_reads_dict[id]["exon3"]))
        exon4MappedReads = str(len(effective_mapped_reads_dict[id]["exon4"]))
        
        countDetail.write(exon2MappedReads + "\t" + exon3MappedReads + "\t" + exon4MappedReads + "\n")
        
        exon2MappedBases = 0
        exon3MappedBases = 0
        exon4MappedBases = 0
        
        for read in effective_mapped_reads_dict[id]["exon2"]:
            exon2MappedBases += effective_mapped_reads_dict[id]["exon2"][read][1]
        for read in effective_mapped_reads_dict[id]["exon3"]:
            exon3MappedBases += effective_mapped_reads_dict[id]["exon3"][read][1]
        for read in effective_mapped_reads_dict[id]["exon4"]:
            exon4MappedBases += effective_mapped_reads_dict[id]["exon4"][read][1]
        
        baseDetail.write(str(exon2MappedBases) + "\t" + str(exon3MappedBases) + "\t" + str(exon4MappedBases) + "\n")
        
    countDetail.close()
    baseDetail.close()
    
    if fileIsNotEmpty(project_folder + "/detail/mappedReadsPerExon.txt") and fileIsNotEmpty(project_folder + "/detail/mappedBasePerExon.txt"):
        return True
    else:
        return False
        
        
################################################################
#
#  Main Function
#
################################################################   
    
def main():
    
    ################################################################
    #
    # Read input parameters
    #
    ################################################################
 
    allele_id_list        = readIdFile(sys.argv[1])
    allele_sequence_list  = readSeqFile(sys.argv[2])
    exons_boundary_dict   = readExonBoundary(sys.argv[3])
    
    project_folder        = sys.argv[4]
    scratch_folder        = sys.argv[5]
    primer_start          = int(sys.argv[6])
    primer_end            = int(sys.argv[7])
    equipment_id          = sys.argv[8]
    
    alignmentLengthCutoff = 50
    
    ################################################################
    #
    # Generate dict for storing effective mapped reads information
    #
    ################################################################
    
    exonWholeCoverage_dict               = {}
    effective_mapped_reads_dict          = {}
    filtered_effective_mapped_reads_dict = {}
    bothPairReadsMapped_dict             = {}
    depthOfCoverage_dict                 = {}
    minimumCoverage_dict                 = {}
    selectedAllele_list                  = []
    alleleCombination_list               = []
    SAM_header_dict                      = {}
    
    ################################################################
    #
    # Go over all allele templates
    #
    ################################################################

    for id_index in range(0, len(allele_id_list)):
    
        output_id               = allele_id_list[id_index]
        output_sequence_length  = len(allele_sequence_list[id_index])
            
        #############################################
        ###
        ###  Check SAM File
        ###
        #############################################
        
        SAM_label  = ""
        SAM_header = []
        
        [SAM_label, SAM_header] = checkSAMFile(scratch_folder, output_id, equipment_id)
        
        if (SAM_label):
            print "%s passed the SAM file checkpoint."%output_id
            SAM_header_dict[output_id] = SAM_header
        else:
            print "%s did not pass the SAM file checkpoint."%output_id
            continue
        
      
        #############################################
        ###
        ###  Extracting alignment info
        ###
        #############################################
        
        [effective_mapped_reads_dict[output_id], exonWholeCoverage_dict[output_id], bothPairReadsMapped_dict[output_id]] = generateExonAlignmentInfo(scratch_folder, output_id, exons_boundary_dict, equipment_id, alignmentLengthCutoff)
    
    #############################################
    ###
    ###  Compute depth of coverage
    ###
    #############################################
    
    [depthOfCoverage_dict, filtered_effective_mapped_reads_dict] = calculateDepthOfCoverage(exons_boundary_dict, effective_mapped_reads_dict, bothPairReadsMapped_dict, scratch_folder, SAM_header_dict)
    
    #############################################
    ###
    ###  Print out depth of coverage
    ###
    #############################################
    
    if printOutDepthOfCoverage(depthOfCoverage_dict, project_folder):
        print "Print depth of coverage: Done!"
    else:
        print "Print depth of coverage failed. please check the output files under folder 'detail'. \n"
    
    #############################################
    ###
    ###  Compute minimum depth of coverage per exon
    ###
    #############################################
    
    minimumCoverage_dict = minimumCoverage(depthOfCoverage_dict)
    
    #############################################
    ###
    ###  Print out minimum coverage
    ###
    #############################################
    
    if printOutMinimumCoverage(minimumCoverage_dict, project_folder):
        print "Print minimum coverage: Done!"
    else:
        print "Print minimum coverage failed. please check the output files under folder 'detail'. \n"
        
    #############################################
    ###
    ###  Print out whole exon coverage
    ###
    #############################################
    
    if printOutWholeExonCoverage(exonWholeCoverage_dict, project_folder):
        print "Print whole exon coverage: Done!"
    else:
        print "Print whole exon coverage failed. please check the output files under folder 'detail'. \n"
    
    
    #############################################
    ###
    ###  Select single cDNA refenrence
    ###
    #############################################
    
    selectedAllele_list = selectSingleAllele(minimumCoverage_dict, exonWholeCoverage_dict, project_folder, primer_start, primer_end)
    
    if len(selectedAllele_list):
        print "We have selected %d alleles\n"%len(selectedAllele_list)
    else:
        print "Something was wrong in selecting single cDNA reference."
        
    #############################################
    ###
    ###  Generate BAM files
    ###
    #############################################
    
    if generateBamFiles(selectedAllele_list, project_folder, scratch_folder):
        print "Successfully generated BAM files!\n"
    else:
        print "Warning: something was wrong in generating BAM files!\n"
        
    #############################################
    ###
    ###  Determine pair-wise allele combination
    ###
    #############################################
    
    # alleleCombination_list[0] includes the mapped base info
    # alleleCombination_list[1] includes the mapped read info
    
    alleleCombination_list = choseAlleleCombination(selectedAllele_list, filtered_effective_mapped_reads_dict, primer_start, primer_end)   
    
    if len(alleleCombination_list[0]):
        print "We have selected %d combos\n"%len(alleleCombination_list[0])
    else:
        print "Something wrong in selecting allele combos."
    
    #############################################
    ###
    ###  Print out mapped reads per exon
    ###
    #############################################
    
    if printOutMappabilityPerExon(filtered_effective_mapped_reads_dict, project_folder):
        print "Printing out mappability per exon: Done!\n"
    else:
        print "Something wrong of <Printing out mappability per exon>, please check the output files under folder 'detail'. \n"
    
    #############################################
    ###
    ###  Print out crossing boundary reads
    ###
    #############################################
    
    #printOutcrossingExonBoundary(effective_mapped_reads_dict, project_folder)
    
    #############################################
    ###
    ###  Print out allele Combinations
    ###
    #############################################
    
    if printOutAlleleCombo(alleleCombination_list[0], project_folder):
        print "Printing Allele Combinations: Done!"
    else:
        print "Something wrong of <Printing Allele Combinations>, please check the output files under folder 'detail'. \n"
    
    #############################################
    ###
    ###  Print out number of mapped reads per allele
    ###
    #############################################
    
    if printOutTotalMappedReads(alleleCombination_list[1], project_folder):
        print "Printing total Mapped Reads: Done!"
    else:
        print "Something wrong of <Printing total Mapped Reads>, please check the output files under folder 'detail'. \n"
   
        
if __name__ == "__main__":
    
    main() 
                
                
                
             
                
                
                     
                    
        
    
