/*******************************************************************************
 * GIGWA - Service implementation
 * Copyright (C) 2016 - 2019, <CIRAD> <IRD>
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License, version 3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * See <http://www.gnu.org/licenses/agpl.html> for details about GNU General
 * Public License V3.
 *******************************************************************************/
package fr.cirad.mgdb.service;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.OutputStream;
import java.lang.reflect.InvocationTargetException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.HttpSession;

import org.apache.avro.AvroRemoteException;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.ga4gh.methods.GAException;
import org.ga4gh.methods.ListReferenceBasesRequest;
import org.ga4gh.methods.ListReferenceBasesResponse;
import org.ga4gh.methods.ReferenceMethods;
import org.ga4gh.methods.SearchCallSetsRequest;
import org.ga4gh.methods.SearchCallSetsResponse;
import org.ga4gh.methods.SearchReferenceSetsRequest;
import org.ga4gh.methods.SearchReferenceSetsResponse;
import org.ga4gh.methods.SearchReferencesRequest;
import org.ga4gh.methods.SearchReferencesResponse;
import org.ga4gh.methods.SearchVariantSetsRequest;
import org.ga4gh.methods.SearchVariantSetsResponse;
import org.ga4gh.methods.SearchVariantsRequest;
import org.ga4gh.methods.VariantMethods;
import org.ga4gh.models.AlleleLocation;
import org.ga4gh.models.AnalysisResult;
import org.ga4gh.models.Call;
import org.ga4gh.models.CallSet;
import org.ga4gh.models.HGVSAnnotation;
import org.ga4gh.models.OntologyTerm;
import org.ga4gh.models.Reference;
import org.ga4gh.models.ReferenceSet;
import org.ga4gh.models.TranscriptEffect;
import org.ga4gh.models.Variant;
import org.ga4gh.models.VariantAnnotation;
import org.ga4gh.models.VariantSet;
import org.ga4gh.models.VariantSetMetadata;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.domain.Sort;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;
import org.springframework.security.core.Authentication;
import org.springframework.security.core.context.SecurityContextHolder;
import org.springframework.stereotype.Component;

import com.mongodb.AggregationOptions;
import com.mongodb.AggregationOptions.Builder;
import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.Bytes;
import com.mongodb.Cursor;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.ServerAddress;

import fr.cirad.controller.GigwaMethods;
import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler;
import fr.cirad.mgdb.importing.SequenceImport;
import fr.cirad.mgdb.importing.VcfImport;
import fr.cirad.mgdb.model.mongo.maintypes.CachedCount;
import fr.cirad.mgdb.model.mongo.maintypes.CustomIndividualMetadata;
import fr.cirad.mgdb.model.mongo.maintypes.CustomIndividualMetadata.CustomIndividualMetadataId;
import fr.cirad.mgdb.model.mongo.maintypes.DBVCFHeader;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.Individual;
import fr.cirad.mgdb.model.mongo.maintypes.Sequence;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData.VariantRunDataId;
import fr.cirad.mgdb.model.mongo.subtypes.AbstractVariantData;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.model.GigwaDensityRequest;
import fr.cirad.model.GigwaSearchReferencesRequest;
import fr.cirad.model.GigwaSearchVariantsExportRequest;
import fr.cirad.model.GigwaSearchVariantsRequest;
import fr.cirad.model.GigwaSearchVariantsResponse;
import fr.cirad.model.GigwaVcfFieldPlotRequest;
import fr.cirad.tools.AlphaNumericComparator;
import fr.cirad.tools.AppConfig;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.GenotypingDataQueryBuilder;
import fr.cirad.tools.mongo.MongoTemplateManager;
import fr.cirad.tools.security.base.AbstractTokenManager;
import fr.cirad.utils.Constants;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;

/**
 *
 * Implementation of Gigwa's data-related functionalities
 *
 * GA4GH / Gigwa concept equivalence:
 * GAReferenceSet = Module / Database, GAVariantSet = Project,
 * GAReference = Sequence, GACallSet = Individual, GAVariant = Variant
 *
 * @author adrien petel, guilhem sempere
 */
@Component
public class GigwaGa4ghServiceImpl implements GigwaMethods, VariantMethods, ReferenceMethods {
    
    /**
     * logger
     */
    protected static final Logger LOG = Logger.getLogger(GigwaGa4ghServiceImpl.class);

    /**
     * The Constant SEQLIST_FOLDER.
     */
    static final public String SEQLIST_FOLDER = "selectedSeqs";

    /**
     * The Constant EXPORT_EXPIRATION_DELAY_MILLIS.
     */
    static final private long EXPORT_EXPIRATION_DELAY_MILLIS = 1000 * 60 * 60 * 24 * 2; /* 2 days */

    /**
     * The Constant TMP_OUTPUT_FOLDER.
     */
    static public final String TMP_OUTPUT_FOLDER = "tmpOutput";

    /**
     * The Constant FRONTEND_URL.
     */
    static final public String FRONTEND_URL = "genofilt";

    static final public String MESSAGE_TEMP_RECORDS_NOT_FOUND = "Unable to find temporary records: please SEARCH again!";

    /**
     * number of simultaneous query threads
     */
    static final private int INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS = 10;
    static final private int MINIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS = 5;
    static final private int MAXIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS = 50;

    static final protected HashMap<String, String> annotationField = new HashMap<>();

	private Boolean fAllowDiskUse = null;
	    
    @Autowired AbstractTokenManager tokenManager;

    @Autowired private AppConfig appConfig;

    /**
     * number format instance
     */
    static protected NumberFormat nf = NumberFormat.getInstance();

    static {
        nf.setMaximumFractionDigits(4);

        annotationField.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + 1, "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + 1);
        annotationField.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + 2, "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + 2);
        annotationField.put(VariantRunData.SECTION_ADDITIONAL_INFO + "." + 1, "$" + VariantRunData.SECTION_ADDITIONAL_INFO);
    }
    
	public boolean isAggregationAllowedToUseDisk() {
		if (appConfig==null) { //We are probably in unit test case
			return true;
		}
		if (fAllowDiskUse == null)
			fAllowDiskUse = !Boolean.parseBoolean(appConfig.get("forbidMongoDiskUse"));
		return fAllowDiskUse.booleanValue();
	}

    @Override
    public LinkedHashSet<String> getProjectVariantTypes(String sModule, int projId) {

        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_VARIANT_TYPES);
       	q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject proj = mongoTemplate.findOne(q, GenotypingProject.class);
        return proj.getVariantTypes();
    }

    @Override
    public List<String> listVariantTypesSorted(String sModule, int projId) {
        List<String> variantTypesArray = new ArrayList(getProjectVariantTypes(sModule, projId));
        Collections.sort(variantTypesArray, new AlphaNumericComparator());
        return variantTypesArray;
    }

    @Override
    public Collection<String> listModules() {
        return MongoTemplateManager.getAvailableModules();
    }

    @Override
    public int getProjectPloidyLevel(String sModule, int projId) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_PLOIDY_LEVEL);
       	q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject proj = mongoTemplate.findOne(q, GenotypingProject.class);
        return proj.getPloidyLevel();
    }

    @Override
    public TreeSet<String> searchableAnnotationFields(String sModule, int projId) {
    	/* This may be more efficient by looking at the VCF header instead */
    	TreeSet<String> result = new TreeSet<>();
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query(Criteria.where("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID).is(projId));
        q.limit(3);
        q.fields().include(VariantRunData.FIELDNAME_SAMPLEGENOTYPES);
        Iterator<VariantRunData> it = mongoTemplate.find(q, VariantRunData.class).iterator();
        while (it.hasNext())
        {
        	VariantRunData vrd = it.next();
            Collection<SampleGenotype> spGTs = vrd.getSampleGenotypes().size() > 100 ? new ArrayList(vrd.getSampleGenotypes().values()).subList(0,  100) : vrd.getSampleGenotypes().values();
            for (HashMap<String, Object> annotationMap : spGTs.stream().map(sg -> sg.getAdditionalInfo()).collect(Collectors.toList()))
            	for (String aiKey : annotationMap.keySet())
 	               if (Number.class.isAssignableFrom(annotationMap.get(aiKey).getClass()))
 	                	result.add(aiKey);
        }
        return result;
    }

    @Override
    public TreeSet<String> getProjectEffectAnnotations(String sModule, int projId) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_EFFECT_ANNOTATIONS);
       	q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject proj = mongoTemplate.findOne(q, GenotypingProject.class);
        return proj.getEffectAnnotations();
    }

    @Override
    public Collection<Integer> getDistinctAlleleCounts(String sModule, Integer projId) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        return mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(GenotypingProject.class)).distinct(GenotypingProject.FIELDNAME_ALLELE_COUNTS, projId == null ? null : new BasicDBObject("_id", projId));
    }

    @Override
    public Map<Integer, String> getProjectIdToNameMap(String sModule) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Map<Integer, String> result = new LinkedHashMap<>();
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_NAME);
        for (GenotypingProject proj : mongoTemplate.find(q, GenotypingProject.class)) {
            result.put(proj.getId(), proj.getName());
        }
        return result;
    }

    @Override
    public String getQueryForGenotypePattern(String gtPattern) {
        return GenotypingDataQueryBuilder.getGenotypePatternToQueryMap().get(gtPattern);
    }

	@Override
	public List<String> listIndividualsInAlphaNumericOrder(String sModule, int project) {
		List<String> indArray = new ArrayList(MgdbDao.getProjectIndividuals(sModule, project));
		Collections.sort(indArray, new AlphaNumericComparator());
		return indArray;
	}

    @Override
    public List<String> listSequences(HttpServletRequest request, String sModule, int projId) {
        List<String> result = getProjectSequences(sModule, projId);

        List<String> externallySelectedSequences = getSequenceIDsBeingFilteredOn(request.getSession(), sModule);
        /* first try to use a list that may have been defined on in a different section of the application (although it may not be limited to the given project) */
        if (externallySelectedSequences != null) {
            result = (List<String>) CollectionUtils.intersection(result, externallySelectedSequences);
        }

        if (result != null) {
            Collections.sort(result, new AlphaNumericComparator());
        }
        return result;
    }

    @Override
    public ArrayList<Object> buildVariantDataQuery(GigwaSearchVariantsRequest gsvr, List<String> externallySelectedSeqs) {
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);
        String actualSequenceSelection = gsvr.getReferenceName();
        if (actualSequenceSelection == null || actualSequenceSelection.length() == 0) {
            if (externallySelectedSeqs != null) {
                actualSequenceSelection = StringUtils.join(externallySelectedSeqs, ";");
            }
        }
        List<String> selectedVariantTypes = gsvr.getSelectedVariantTypes().length() == 0 ? null : Arrays.asList(gsvr.getSelectedVariantTypes().split(";"));
        List<String> selectedSequences = Arrays.asList(actualSequenceSelection == null || actualSequenceSelection.length() == 0 ? new String[0] : actualSequenceSelection.split(";"));
        List<String> alleleCountList = gsvr.getAlleleCount().length() == 0 ? null : Arrays.asList(gsvr.getAlleleCount().split(";"));
        
        BasicDBList variantFeatureFilterList = new BasicDBList();

        /* Step to match selected variant types */
        if (selectedVariantTypes != null && selectedVariantTypes.size() > 0) {
            BasicDBList orList1 = new BasicDBList();
            DBObject orSelectedVariantTypesList = new BasicDBObject();
            for (String aSelectedVariantTypes : selectedVariantTypes) {
                DBObject orClause1 = new BasicDBObject(VariantData.FIELDNAME_TYPE, aSelectedVariantTypes);
                orList1.add(orClause1);
                orSelectedVariantTypesList.put("$or", orList1);
            }
            variantFeatureFilterList.add(orSelectedVariantTypesList);
        }

        /* Step to match selected chromosomes */
        if (selectedSequences != null && selectedSequences.size() > 0 && selectedSequences.size() != getProjectSequences(sModule, projId).size()) {
            variantFeatureFilterList.add(new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, new BasicDBObject("$in", selectedSequences)));
        }

        /* Step to match variants that have a position included in the specified range */
        if (gsvr.getStart() != null || gsvr.getEnd() != null) {
            if (gsvr.getStart() != null && gsvr.getStart() != -1) {
                DBObject firstPosStart = new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, new BasicDBObject("$gte", gsvr.getStart()));
                variantFeatureFilterList.add(firstPosStart);
            }
            if (gsvr.getEnd() != null && gsvr.getEnd() != -1) {
                DBObject lastPosStart = new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, new BasicDBObject("$lte", gsvr.getEnd()));
                variantFeatureFilterList.add(lastPosStart);
            }
        }

        /* Step to match selected number of alleles */
        if (alleleCountList != null) {
            BasicDBList orList3 = new BasicDBList();
            DBObject orSelectedNumberOfAllelesList = new BasicDBObject();
            for (String aSelectedNumberOfAlleles : alleleCountList) {
                int alleleNumber = Integer.parseInt(aSelectedNumberOfAlleles);
                orList3.add(new BasicDBObject(VariantData.FIELDNAME_KNOWN_ALLELE_LIST, new BasicDBObject("$size", alleleNumber)));
                orSelectedNumberOfAllelesList.put("$or", orList3);
            }
            variantFeatureFilterList.add(orSelectedNumberOfAllelesList);
        }

        return variantFeatureFilterList;
    }
    
    private String isSearchedDatasetReasonablySized(GigwaSearchVariantsRequest gsvr)
    {
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);
        
        int nProjectIndCount = MgdbDao.getProjectIndividuals(sModule, projId).size();
        int nGroup1IndCount = gsvr.getCallSetIds() != null && gsvr.getCallSetIds().size() != 0 ? gsvr.getCallSetIds().size() : nProjectIndCount;
        int nGroup2IndCount = gsvr.getCallSetIds2() != null && gsvr.getCallSetIds2().size() != 0 ? gsvr.getCallSetIds().size() : nProjectIndCount;
        
    	List<Integer> groupsForWhichToFilterOnGenotypingData = GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingData(gsvr);
    	int nIndCount = (groupsForWhichToFilterOnGenotypingData.contains(0) ? nGroup1IndCount : 0)  + (groupsForWhichToFilterOnGenotypingData.contains(1) ? nGroup2IndCount : 0);
    	if (nIndCount == 0)
    		return null;	// no genotyping data filtering involved

    	int nMaxBillionGenotypesInvolved = 1;	// default
    	try
    	{
    		nMaxBillionGenotypesInvolved = Integer.parseInt(appConfig.get("maxSearchableBillionGenotypes"));
    	}
    	catch (Exception ignored)
    	{}

        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_SEQUENCES);
       	q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject proj = mongoTemplate.findOne(q, GenotypingProject.class);
        
    	int nSelectedSeqCount = gsvr.getReferenceName() == null || gsvr.getReferenceName().length() == 0 ? proj.getSequences().size() : gsvr.getReferenceName().split(";").length;
    	if (nSelectedSeqCount == 1)
    		return null;	// we can't expect user to select less than a single sequence

    	int nAvgVariantsPerSeq = (int) (mongoTemplate.count(null, VariantData.class) / Math.max(1, proj.getSequences().size()));
    	BigInteger maxSeqCount = BigInteger.valueOf(1000000000).multiply(BigInteger.valueOf(nMaxBillionGenotypesInvolved)).divide(BigInteger.valueOf(nAvgVariantsPerSeq).multiply(BigInteger.valueOf(nIndCount)));
    	int nMaxSeqCount = Math.max(1, maxSeqCount.intValue());

    	if (nSelectedSeqCount <= nMaxSeqCount)
    		return null;
    	
    	return "This is a big database. Given the number of selected individuals, you may only work on a maximum of " + nMaxSeqCount + " sequence(s) at a time.";
    }

    @Override
    public long countVariants(GigwaSearchVariantsRequest gsvr, boolean fSelectionAlreadyExists) throws Exception {
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        boolean fGotTokenManager = tokenManager != null;	// if null, we are probably being invoked via unit-test
        String token = !fGotTokenManager ? Helper.convertToMD5(String.valueOf(System.currentTimeMillis())) /* create a mock up token */ : tokenManager.readToken(gsvr.getRequest());

        ProgressIndicator progress = ProgressIndicator.get(token);
        if (progress == null) {
            progress = new ProgressIndicator(token, new String[0]);
            ProgressIndicator.registerProgressIndicator(progress);
        }
    	String sizeProblemMsg = gsvr.shallApplyMatrixSizeLimit() ? isSearchedDatasetReasonablySized(gsvr) : null;
    	if (sizeProblemMsg != null)
    	{
            progress.setError(sizeProblemMsg);
            return 0;
    	}

        DBCollection tmpVarColl = getTemporaryVariantCollection(sModule, progress.getProcessId(), fGotTokenManager && !fSelectionAlreadyExists);
        String queryKey = getQueryKey(gsvr);

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        DBCollection cachedCountCollection = mongoTemplate.getCollection(mongoTemplate.getCollectionName(CachedCount.class));
//			cachedCountCollection.drop();
        DBCursor countCursor = cachedCountCollection.find(new BasicDBObject("_id", queryKey));
        Long count = null;
        if (countCursor.hasNext()) {
            count = 0l;
            for (Object aPartialCount : ((BasicDBList) countCursor.next().get(MgdbDao.FIELD_NAME_CACHED_COUNT_VALUE)).toArray()) {
                count += (Long) aPartialCount;
            }
        }
        LOG.debug((count == null ? "new" : "existing") + " queryKey hash: " + queryKey);
        if (count == null)
        {
            long before = System.currentTimeMillis();
            progress.addStep("Counting matching variants");
            
            List<String> alleleCountList = gsvr.getAlleleCount().length() == 0 ? null : Arrays.asList(gsvr.getAlleleCount().split(";"));

            GenotypingProject genotypingProject = mongoTemplate.findById(projId, GenotypingProject.class);
            if (genotypingProject.getAlleleCounts().size() != 1 || genotypingProject.getAlleleCounts().iterator().next() != 2) {	// Project does not only have bi-allelic data: make sure we can apply MAF filter on selection
                boolean fExactlyOneNumberOfAllelesSelected = alleleCountList != null && alleleCountList.size() == 1;
                boolean fBiAllelicSelected = fExactlyOneNumberOfAllelesSelected && "2".equals(alleleCountList.get(0));
                boolean fMafRequested = (gsvr.getMaxmaf() != null && gsvr.getMaxmaf() < 50) || (gsvr.getMinmaf() != null && gsvr.getMinmaf() > 0);
                if (fMafRequested && !fBiAllelicSelected) {
                    progress.setError("MAF is only supported on biallelic data!");
                    countCursor.close();
                    return 0l;
                }
            }

            DBCollection varColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
            List<Integer> filteredGroups = GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsvr);
            BasicDBList variantQueryDBList = (BasicDBList) buildVariantDataQuery(gsvr, !fGotTokenManager ? null : getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), sModule));
                            
            if (variantQueryDBList.isEmpty()) {
                if (filteredGroups.size() == 0 && mongoTemplate.count(null, GenotypingProject.class) == 1)
                    count = mongoTemplate.count(new Query(), VariantData.class);	// no filter whatsoever
            }
            else if (filteredGroups.size() == 0) {	// filtering on variant features only: we just need a count
                count = varColl.count(new BasicDBObject("$and", variantQueryDBList));
            }

            if (count != null) {
                BasicDBObject dbo = new BasicDBObject("_id", queryKey);
                dbo.append(MgdbDao.FIELD_NAME_CACHED_COUNT_VALUE, new Long[]{count});
                cachedCountCollection.save(dbo);
            }
            else
            {	// filter on genotyping data
                boolean fPreFilterOnVarColl = false;
                List<ServerAddress> mongoServerList = mongoTemplate.getDb().getMongo().getServerAddressList();
                boolean fMongoOnSameServer = mongoServerList.size() == 1 && Arrays.asList("127.0.0.1", "localhost").contains(mongoServerList.get(0).getHost());
                if (variantQueryDBList.size() > 0)
                {
                	Number avgObjSize = (Number) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
                    if (avgObjSize.doubleValue() >= 10240)
                    {	// it may be worth pre-filtering data on variant collection because filtering speed on the run collection is affected by the document size
    	                long totalCount = mongoTemplate.count(new Query(), VariantData.class), preFilterCount = varColl.count(new BasicDBObject("$and", variantQueryDBList));
    	                fPreFilterOnVarColl = preFilterCount <= totalCount*(fMongoOnSameServer ? .85 : .45);	// only pre-filter if less than a given portion of the total variants are to be retained
    	                if (fPreFilterOnVarColl)
    	                	LOG.debug("Pre-filtering data on variant collection");
                    }
                }
                GenotypingDataQueryBuilder genotypingDataQueryBuilder = new GenotypingDataQueryBuilder(gsvr, tmpVarColl, variantQueryDBList, true);
            	final int nChunkCount = genotypingDataQueryBuilder.getNumberOfQueries();
            	final List<Integer> shuffledChunkIndexes = genotypingDataQueryBuilder.suffleChunkOrder();

                try
                {
                    if (nChunkCount > 1)
                        LOG.debug("Query split into " + nChunkCount);

                    final Long[] partialCountArray = new Long[nChunkCount];
                    final Builder aggOpts = AggregationOptions.builder().allowDiskUse(isAggregationAllowedToUseDisk());
                    final ArrayList<Thread> threadsToWaitFor = new ArrayList<>();
                    final AtomicInteger finishedThreadCount = new AtomicInteger(0);

                    int i = -1, nNConcurrentThreads = INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS;
                    while (genotypingDataQueryBuilder.hasNext()) {
                        final List<DBObject> genotypingDataPipeline = genotypingDataQueryBuilder.next();

                        final int chunkIndex = shuffledChunkIndexes.get(++i);
                        
                        boolean fMultiProjectDB = false;

                        BasicDBObject initialMatch = (BasicDBObject) genotypingDataPipeline.get(0).get("$match");
                        if (initialMatch != null && fPreFilterOnVarColl)
                        {	// modifiedMatchAnd will be the one applied to variants collection when pre-filtering
                        	BasicDBList modifiedMatchAnd = (BasicDBList) ((BasicDBList) initialMatch.get("$and")).clone();
                        	if (modifiedMatchAnd != null)
                        	{
                        		List<DBObject> toAdd = new ArrayList<>(), toRemove = new ArrayList<>();
                        		for (Object filter : modifiedMatchAnd)
                        		{
                        			Object variantIdFilter = ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
                        			if (variantIdFilter != null)
                        			{
                        				toAdd.add(new BasicDBObject("_id", variantIdFilter));
                        				toRemove.add((DBObject) filter);
                        			}
                        			else if (null != ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID)) {
                        				toRemove.add((DBObject) filter);	// no project info to filter on in the variants collection
	                    				fMultiProjectDB = true;
	                    			}

                        		}
                        		modifiedMatchAnd.addAll(toAdd);
                        		modifiedMatchAnd.removeAll(toRemove);
                        	}

                        	if (fMongoOnSameServer)
                        	{	// always worth pre-filtering
                                DBCursor variantCursor = varColl.find(new BasicDBObject("$and", modifiedMatchAnd), new BasicDBObject("_id", 1));
                                List<Comparable> chunkPreFilteredIDs = new ArrayList<>();
                                while (variantCursor.hasNext())
                                	chunkPreFilteredIDs.add((Comparable) variantCursor.next().get("_id"));
                                if (chunkPreFilteredIDs.size() == 0)
                                {	// no variants match indexed part of the query: skip chunk
                                	partialCountArray[chunkIndex] = 0l;
                                	// do as if one more async thread was launched so we keep better track of the progress
                                	threadsToWaitFor.add(null);
                                	finishedThreadCount.incrementAndGet();
                                	continue;
                                }
                                else
                                {	// DB server is the same machine as web server: $in operator will not be expensive
                                	if (!fMultiProjectDB)	// for single project dbs, $in is equivalent to original query, otherwise only a pre-filter
                                		genotypingDataPipeline.remove(0);
        	                        genotypingDataPipeline.add(0, new BasicDBObject("$match", new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", chunkPreFilteredIDs))));
                                }
                        	}
                        	else
                        	{	// only try and use pre-filtering to avoid executing genotyping data queries on irrelevant chunks
                                if (varColl.count(new BasicDBObject("$and", modifiedMatchAnd)) == 0)
                                {	// no variants match indexed part of the query: skip chunk
                                	partialCountArray[chunkIndex] = 0l;
                                	// do as if one more async thread was launched so we keep better track of the progress
                                	threadsToWaitFor.add(null);
                                	finishedThreadCount.incrementAndGet();
                                	continue;
                                }
                        	}
                        }

                        // Now the $group operation, used for counting
                        genotypingDataPipeline.add(new BasicDBObject("$count", "count"));

                        if (progress.isAborted()) {
                            return 0l;
                        }

                        final ProgressIndicator finalProgress = progress;

                        Thread queryThread = new Thread() {
                            @Override
                            public void run() {
                                Cursor it = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(genotypingDataPipeline, aggOpts.build());
                                partialCountArray[chunkIndex] = it.hasNext() ? ((Number) it.next().get("count")).longValue() : 0;
                                finalProgress.setCurrentStepProgress((short) (finishedThreadCount.incrementAndGet() * 100 / nChunkCount));
                                genotypingDataPipeline.clear();	// release memory (VERY IMPORTANT)
                                it.close();
                            }
                        };

                        if (chunkIndex % nNConcurrentThreads == (nNConcurrentThreads - 1)) {
                            threadsToWaitFor.add(queryThread); // only needed to have an accurate count
                            queryThread.run();	// run synchronously
                            
                            // regulate number of concurrent threads
                            int nRunningThreadCount = threadsToWaitFor.size() - finishedThreadCount.get();
                            if (nRunningThreadCount > nNConcurrentThreads * .5)
                            	nNConcurrentThreads = (int) (nNConcurrentThreads / 1.5);
                            else if (nRunningThreadCount < nNConcurrentThreads * .25)
                            	nNConcurrentThreads *= 1.5;
                            nNConcurrentThreads = Math.min(MAXIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS, Math.max(MINIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS, nNConcurrentThreads));
                        }
                        else {
                            threadsToWaitFor.add(queryThread);
                            queryThread.start();	// run asynchronously for better speed
                        }
                    }

                    for (Thread t : threadsToWaitFor) // wait for all threads before moving to next phase
                    	if (t != null)
                    		t.join();
                    progress.setCurrentStepProgress(100);
                    
                    count = 0l;
                    for (Long partialCount : partialCountArray) {
                        count += partialCount;
                    }

                    BasicDBObject dbo = new BasicDBObject("_id", queryKey);
                    dbo.append(MgdbDao.FIELD_NAME_CACHED_COUNT_VALUE, partialCountArray);
                    cachedCountCollection.save(dbo);
                }
                catch (InterruptedException e) {
                    LOG.debug("InterruptedException", e);
                }
            }
            LOG.info("countVariants found " + count + " results in " + (System.currentTimeMillis() - before) / 1000d + "s");
        }
        
	if (!fSelectionAlreadyExists)
	    progress.markAsComplete();
        if (progress.isAborted()) {
            return 0l;
        }
        return count;
    }

	/**
     * Gets the temporary variant collection.
     *
     * @param sModule the module
     * @param processID the process id
     * @param fEmptyItBeforeHand whether or not to empty it beforehand
     * @return the temporary variant collection
     */
    private DBCollection getTemporaryVariantCollection(String sModule, String processID, boolean fEmptyItBeforeHand) {
        DBCollection tmpColl = MongoTemplateManager.get(sModule).getCollection(MongoTemplateManager.TEMP_COLL_PREFIX + Helper.convertToMD5(processID));
        if (fEmptyItBeforeHand) {
        	
//        	ArrayList<StackTraceElement> keptStackTraceElements = new ArrayList<>();
//        	Exception e = new Exception("Check stack trace");
//        	for (StackTraceElement ste : e.getStackTrace())
//        		if (ste.toString().startsWith("fr.cirad."))
//        			keptStackTraceElements.add(ste);
//        	e.setStackTrace(keptStackTraceElements.toArray(new StackTraceElement[keptStackTraceElements.size()]));
//            LOG.debug("Dropping " + sModule + "." + tmpColl.getName() + " from getTemporaryVariantCollection", e);
            
            tmpColl.drop();
    		BasicDBObject runCollIndexKeys = new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, 1);
    		runCollIndexKeys.put(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, 1);
    		tmpColl.createIndex(runCollIndexKeys);
            tmpColl.createIndex(VariantData.FIELDNAME_TYPE);
        }
        return tmpColl;
    }

    /**
     * get a cursor of variant records and get a given number of result from it
     * @throws Exception 
     */
    private DBCursor getVariantCursor(GigwaSearchVariantsRequest gsvr) throws Exception {

        DBCursor cursor = null;
        String token = tokenManager.readToken(gsvr.getRequest());
        
        if (token == null)
            return null;
        
        ProgressIndicator progress = ProgressIndicator.get(token);
        if (progress == null) {
            progress = new ProgressIndicator(token, new String[0]);
            ProgressIndicator.registerProgressIndicator(progress);
        }

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];

        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);

        DBCollection variantColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
        DBCollection tempVarColl = getTemporaryVariantCollection(sModule, token, false);
        BasicDBList variantQueryDBList = (BasicDBList) buildVariantDataQuery(gsvr, getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), sModule));

        DBCollection varCollForBuildingRows = tempVarColl.count() == 0 ? variantColl : tempVarColl;
        cursor = varCollForBuildingRows.find(!variantQueryDBList.isEmpty() ? new BasicDBObject("$and", variantQueryDBList) : null);
        if (gsvr.getSortBy() != null && gsvr.getSortBy().length() > 0)
            cursor.sort(new BasicDBObject(gsvr.getSortBy(), Integer.valueOf("DESC".equalsIgnoreCase(gsvr.getSortDir()) ? -1 : 1)));
        else
        {
        	BasicDBObject sortObj = new BasicDBObject(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, 1);
        	sortObj.put(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, 1);
            cursor.sort(sortObj);
        }
        // skip the results we don't want           
        cursor.skip(Integer.parseInt(gsvr.getPageToken()) * gsvr.getPageSize()).limit(gsvr.getPageSize());

        progress.markAsComplete();
        return cursor;
    }

    @Override
    public long findVariants(GigwaSearchVariantsRequest gsvr) throws Exception {
        String token = tokenManager.readToken(gsvr.getRequest());

        final ProgressIndicator progress = new ProgressIndicator(token, new String[0]);
        ProgressIndicator.registerProgressIndicator(progress);
        String sizeProblemMsg = gsvr.shallApplyMatrixSizeLimit() ? isSearchedDatasetReasonablySized(gsvr) : null;
    	if (sizeProblemMsg != null)
    	{
            progress.setError(sizeProblemMsg);
            return 0;
    	}
    	
        progress.addStep("Finding matching variants");

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        String queryKey = getQueryKey(gsvr);
        final DBCollection tmpVarColl = getTemporaryVariantCollection(sModule, progress.getProcessId(), true);
        
        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        DBCollection cachedCountCollection = mongoTemplate.getCollection(mongoTemplate.getCollectionName(CachedCount.class));
        DBCursor countCursor = cachedCountCollection.find(new BasicDBObject("_id", queryKey));

        final Object[] partialCountArray = !countCursor.hasNext() ? null : ((BasicDBList) countCursor.next().get(MgdbDao.FIELD_NAME_CACHED_COUNT_VALUE)).toArray();
        final LinkedHashMap<Integer, Long> partialCountMap = new LinkedHashMap<>();	// progress display will be more accurate if we skip empty chunks
        long nTotalCount = 0;
        if (partialCountArray != null)
        {
        	for (int i=0; i<partialCountArray.length; i++)
        		if ((long) partialCountArray[i] != 0)
        		{
        			long n = (long) partialCountArray[i];
        			partialCountMap.put(i, n);
        			nTotalCount += n;
        		}
	        if (nTotalCount == 0)
	        {
	            progress.markAsComplete();
	        	return 0;
	        }
        }

        List<Integer> filteredGroups = GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsvr);
        BasicDBList variantQueryDBList = (BasicDBList) buildVariantDataQuery(gsvr, getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), sModule));

    	long before = System.currentTimeMillis();
        if (filteredGroups.size() > 0)
        {	// filter on genotyping data
        	final ArrayList<Thread> threadsToWaitFor = new ArrayList<>();
            final AtomicInteger finishedThreadCount = new AtomicInteger(0);
            final GenotypingDataQueryBuilder genotypingDataQueryBuilder = new GenotypingDataQueryBuilder(gsvr, tmpVarColl, variantQueryDBList, false);
            try
            {
            	final int nChunkCount = genotypingDataQueryBuilder.getNumberOfQueries();
            	final List<Integer> shuffledChunkIndexes = genotypingDataQueryBuilder.suffleChunkOrder();
            	
            	final Long[] partialCountArrayToFill = partialCountArray == null ? new Long[nChunkCount] : null;
            	if (partialCountArrayToFill != null)
            		LOG.info("Find without prior count: do both at once");
            	else if (nChunkCount != partialCountArray.length) {
                    progress.setError("Different number of chunks between counting and listing variant rows!");
                    LOG.error(progress.getError());
                    return 0;
                }
            	
                DBCollection varColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
                
                boolean fPreFilterOnVarColl = false;
                List<ServerAddress> mongoServerList = mongoTemplate.getDb().getMongo().getServerAddressList();
                boolean fMongoOnSameServer = mongoServerList.size() == 1 && Arrays.asList("127.0.0.1", "localhost").contains(mongoServerList.get(0).getHost());
                if (variantQueryDBList.size() > 0)
                {
                	Number avgObjSize = (Number) mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).getStats().get("avgObjSize");
                    if (avgObjSize.doubleValue() >= 10240)
                    {	// it may be worth pre-filtering data on variant collection because filtering speed on the run collection is affected by the document size
    	                long totalCount = mongoTemplate.count(new Query(), VariantData.class), preFilterCount = varColl.count(new BasicDBObject("$and", variantQueryDBList));
    	                fPreFilterOnVarColl = preFilterCount <= totalCount*(fMongoOnSameServer ? .85 : .45);	// only pre-filter if less than a given portion of the total variants are to be retained
    	                if (fPreFilterOnVarColl)
    	                	LOG.debug("Pre-filtering data on variant collection");
                    }
                }

                if (nChunkCount > 1)
                    LOG.debug("Query split into " + nChunkCount);

                int i = -1, nNConcurrentThreads = INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS;
                while (genotypingDataQueryBuilder.hasNext()) {
                    final List<DBObject> genotypingDataPipeline = genotypingDataQueryBuilder.next();
                    if (progress.isAborted() || progress.getError() != null)
                        return 0;

                    final int chunkIndex = shuffledChunkIndexes.get(++i);                    
                    if (partialCountMap.size() > 0 && !partialCountMap.containsKey(chunkIndex))
                        continue;
                    
                    boolean fMultiProjectDB = false;

                    BasicDBObject initialMatch = (BasicDBObject) genotypingDataPipeline.get(0).get("$match");
                    if (initialMatch != null && fPreFilterOnVarColl)
                    {
                    	BasicDBList modifiedMatchAnd = (BasicDBList) ((BasicDBList) initialMatch.get("$and")).clone();
                    	if (modifiedMatchAnd != null)
                    	{
                    		List<DBObject> toAdd = new ArrayList<>(), toRemove = new ArrayList<>();
                    		for (Object filter : modifiedMatchAnd)
                    		{
                    			Object variantIdFilter = ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
                    			if (variantIdFilter != null)
                    			{
                    				toAdd.add(new BasicDBObject("_id", variantIdFilter));
                    				toRemove.add((DBObject) filter);
                    			}
                    			else if (null != ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID)) {
                    				toRemove.add((DBObject) filter);	// no project info to filter on in the variants collection
                    				fMultiProjectDB = true;
                    			}
                    		}
                    		modifiedMatchAnd.addAll(toAdd);
                    		modifiedMatchAnd.removeAll(toRemove);
                    	}

                    	if (fMongoOnSameServer)
                    	{	// always worth pre-filtering
                            DBCursor variantCursor = varColl.find(new BasicDBObject("$and", modifiedMatchAnd), new BasicDBObject("_id", 1));
                            List<Comparable> chunkPreFilteredIDs = new ArrayList<>();
                            while (variantCursor.hasNext())
                            	chunkPreFilteredIDs.add((Comparable) variantCursor.next().get("_id"));
                            if (chunkPreFilteredIDs.size() == 0)
                            {	// no variants match indexed part of the query: skip chunk
	                            if (partialCountArrayToFill != null)
	                            	partialCountArrayToFill[chunkIndex] = 0l;
                            	// do as if one more async thread was launched so we keep better track of the progress
                            	threadsToWaitFor.add(null);
                            	finishedThreadCount.incrementAndGet();
                            	continue;
                            }
                            else
                            {	// DB server is the same machine as web server: $in operator will not be expensive
                            	if (!fMultiProjectDB)	// for single project dbs, $in is equivalent to original query, otherwise only a pre-filter
                            		genotypingDataPipeline.remove(0);
    	                        genotypingDataPipeline.add(0, new BasicDBObject("$match", new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", chunkPreFilteredIDs))));
                            }
                    	}
                    	else
                    	{	// only try and use pre-filtering to avoid executing genotyping data queries on irrelevant chunks
                            if (varColl.count(new BasicDBObject("$and", modifiedMatchAnd)) == 0)
                            {	// no variants match indexed part of the query: skip chunk
	                            if (partialCountArrayToFill != null)
	                            	partialCountArrayToFill[chunkIndex] = 0l;
                            	// do as if one more async thread was launched so we keep better track of the progress
                            	threadsToWaitFor.add(null);
                            	finishedThreadCount.incrementAndGet();
                            	continue;
                            }
                    	}
                    }

                    Thread queryThread = new Thread() {
                        @Override
                        public void run() {
                        	try {
                        		if (partialCountArray != null)
                        			genotypingDataPipeline.add(new BasicDBObject("$limit", partialCountArray[chunkIndex]));
                        		
	                            Cursor genotypingDataCursor = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(genotypingDataPipeline, AggregationOptions.builder().allowDiskUse(isAggregationAllowedToUseDisk()).build());
	                            final ArrayList<DBObject> variantsThatPassedRunFilterForThisChunk = new ArrayList<>();
	                            while (genotypingDataCursor.hasNext())
	                            {
	                            	DBObject variant = genotypingDataCursor.next();
	                            	variant.removeField(GenotypingDataQueryBuilder.MAIN_RESULT_PROJECTION_FIELD);
	                            	variant.removeField(GenotypingDataQueryBuilder.MOSTLY_SAME_GENOTYPE_RESULT_PROJECTION_FIELD);
	                                variantsThatPassedRunFilterForThisChunk.add(variant);
	                            }
	
	                            if (partialCountArrayToFill != null)
	                            	partialCountArrayToFill[chunkIndex] = (long) variantsThatPassedRunFilterForThisChunk.size();
	                            if (variantsThatPassedRunFilterForThisChunk.size() > 0)
	                            	tmpVarColl.insert(variantsThatPassedRunFilterForThisChunk);
	                            
	                            genotypingDataPipeline.clear();	// release memory (VERY IMPORTANT)
	                            genotypingDataCursor.close();
                        	}
                        	catch (Exception e) {
                        		LOG.error("Error", e);
                        		progress.setError(e.getMessage());
                        	}
                        }
                    };

                    if (chunkIndex % nNConcurrentThreads == (nNConcurrentThreads - 1)) {
                        threadsToWaitFor.add(queryThread); // we only need to have an accurate count
                        queryThread.run();	// run synchronously
                        
                        // regulate number of concurrent threads
                        int nRunningThreadCount = threadsToWaitFor.size() - finishedThreadCount.get();
                        if (nRunningThreadCount > nNConcurrentThreads * .5)
                        	nNConcurrentThreads = (int) (nNConcurrentThreads / 1.5);
                        else if (nRunningThreadCount < nNConcurrentThreads * .25)
                        	nNConcurrentThreads *= 1.5;
                        nNConcurrentThreads = Math.min(MAXIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS, Math.max(MINIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS, nNConcurrentThreads));
//                        System.out.println(nRunningThreadCount + " / " + threadsToWaitFor.size() + " -> " + nNConcurrentThreads);
                    }
                    else {
                        threadsToWaitFor.add(queryThread);
                        queryThread.start();	// run asynchronously for better speed
                    }
                    progress.setCurrentStepProgress((short) (i * 100 / nChunkCount));
                }

                for (Thread t : threadsToWaitFor) // wait for all threads before moving to next phase
                	if (t != null)
                		t.join();
                progress.setCurrentStepProgress(100);
                
                if (partialCountArrayToFill != null)
                {
	                BasicDBObject dbo = new BasicDBObject("_id", queryKey);
	                dbo.append(MgdbDao.FIELD_NAME_CACHED_COUNT_VALUE, partialCountArrayToFill);
	                cachedCountCollection.save(dbo);
                }
            }
            catch (InterruptedException e) {
                LOG.debug("InterruptedException : " + e);
                // throw e;
            }
        }
        
        if (partialCountArray == null)
        	nTotalCount = countVariants(gsvr, true);

        if (progress.isAborted() || progress.getError() != null)
            return 0;

        progress.markAsComplete();
        LOG.info("findVariants found " + nTotalCount + " results in " + (System.currentTimeMillis() - before) / 1000d + "s");
        return nTotalCount;
    }

    @Override
    public void exportVariants(GigwaSearchVariantsExportRequest gsver, String token, HttpServletResponse response) throws Exception
    {
    	String processId = "export_" + token;
        final ProgressIndicator progress = new ProgressIndicator(processId, new String[0]); 
        ProgressIndicator.registerProgressIndicator(progress);
        
        new Thread() { public void run() {
        		try {
					cleanupExpiredExportData(gsver.getRequest());
				}
        		catch (IOException e) {
					LOG.error("Unable to cleanup expired export files", e);
				}
        	}
        }.start();

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsver.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        long before = System.currentTimeMillis();
        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
    	int nGroupsToFilterGenotypingDataOn = GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsver).size();

		List<String> selectedIndividualList1 = gsver.getCallSetIds().size() == 0 ? MgdbDao.getProjectIndividuals(sModule, projId) /* no selection means all selected */ : gsver.getCallSetIds();
		List<String> selectedIndividualList2 = gsver.getCallSetIds2().size() == 0 && nGroupsToFilterGenotypingDataOn > 1 ? MgdbDao.getProjectIndividuals(sModule, projId) /* no selection means all selected */ : gsver.getCallSetIds2();
		Collection<String> individualsToExport = gsver.getExportedIndividuals();
		if (individualsToExport.size() == 0)
			individualsToExport = MgdbDao.getProjectIndividuals(sModule, projId);

        long count = countVariants(gsver, true);
        DBCollection tmpVarColl = getTemporaryVariantCollection(sModule, token, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getName());
        final BasicDBList variantQueryDBList = (BasicDBList) buildVariantDataQuery(gsver, getSequenceIDsBeingFilteredOn(gsver.getRequest().getSession(), sModule));

		if (nGroupsToFilterGenotypingDataOn > 0 && nTempVarCount == 0)
		{
			progress.setError(MESSAGE_TEMP_RECORDS_NOT_FOUND);
			return;
		}

        // use a cursor to avoid using too much memory
        String sequenceField = VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE;
        String startField = VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE;
    	BasicDBObject sortObj = new BasicDBObject(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, 1);
    	sortObj.put(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, 1);
        DBObject projection = new BasicDBObject();
        projection.put(sequenceField, 1);
        projection.put(startField, 1);

        String mainVarCollName = mongoTemplate.getCollectionName(VariantData.class), usedVarCollName = nTempVarCount == 0 ? mainVarCollName : tmpVarColl.getName();
        DBCollection usedVarColl = mongoTemplate.getCollection(usedVarCollName);
        BasicDBObject variantQuery = nTempVarCount == 0 && !variantQueryDBList.isEmpty() ? new BasicDBObject("$and", variantQueryDBList) : null;
        if (gsver.shallApplyMatrixSizeLimit())
    	{	// make sure the matrix is not too big
        	int nMaxBillionGenotypesInvolved = 1;	// default
        	try {
                Authentication auth = tokenManager.getAuthenticationFromToken(tokenManager.readToken(gsver.getRequest()));
        		nMaxBillionGenotypesInvolved = Integer.parseInt(appConfig.get("maxExportableBillionGenotypes_" + (auth == null ? "anonymousUser" : auth.getName())));
        	}
        	catch (Exception ignored1) {
               	try {
            		nMaxBillionGenotypesInvolved = Integer.parseInt(appConfig.get("maxExportableBillionGenotypes"));
            	}
            	catch (Exception ignored2)
            	{}
        	}
        	
        	if (nMaxBillionGenotypesInvolved == 0)
        	{
	            progress.setError("You are not allowed to export any genotyping data.");
	            return;
        	}
        	
        	BigInteger matrixSize = BigInteger.valueOf(usedVarColl.count(variantQuery)).multiply(BigInteger.valueOf(individualsToExport.size()));
        	BigInteger maxAllowedSize = BigInteger.valueOf(1000000000).multiply(BigInteger.valueOf(nMaxBillionGenotypesInvolved));
        	
        	if (matrixSize.divide(maxAllowedSize).intValue() >= 1)
        	{
	            progress.setError("You may only export up to " + nMaxBillionGenotypesInvolved + " billion genotypes. The current selection contains " + BigDecimal.valueOf(matrixSize.longValue()).divide(BigDecimal.valueOf(1000000000)).setScale(2, BigDecimal.ROUND_HALF_UP) + " billion genotypes.");
	            return;
        	}
    	}

    	progress.addStep("Identifying matching variants");
        DBCursor markerCursor = usedVarColl.find(variantQuery, projection).sort(sortObj);
        markerCursor.addOption(Bytes.QUERYOPTION_NOTIMEOUT);
        OutputStream os = null;

        try
        {
			AbstractIndividualOrientedExportHandler individualOrientedExportHandler = AbstractIndividualOrientedExportHandler.getIndividualOrientedExportHandlers().get(gsver.getExportFormat());
			AbstractMarkerOrientedExportHandler markerOrientedExportHandler = AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().get(gsver.getExportFormat());

            String filename = sModule + "__project" + projId + "__" + new SimpleDateFormat("yyyy-MM-dd").format(new Date()) + "__" + count + "variants__" + gsver.getExportFormat() + "." + (individualOrientedExportHandler != null ? individualOrientedExportHandler : markerOrientedExportHandler).getExportArchiveExtension();

            LOG.info((gsver.isKeepExportOnServer() ? "On-server" : "Direct-download") + " export requested: " + processId);
            if (gsver.isKeepExportOnServer()) {
            	Authentication auth = SecurityContextHolder.getContext().getAuthentication();
            	String sExportingUser = auth == null || "anonymousUser".equals(auth.getName()) ? "anonymousUser" : auth.getName();
                String relativeOutputFolder = FRONTEND_URL + File.separator + TMP_OUTPUT_FOLDER + File.separator + sExportingUser + File.separator + Helper.convertToMD5(processId) + File.separator;
                File outputLocation = new File(gsver.getRequest().getSession().getServletContext().getRealPath(File.separator + relativeOutputFolder));
                if (!outputLocation.exists() && !outputLocation.mkdirs()) {
                    throw new Exception("Unable to create folder: " + outputLocation);
                }
                os = new FileOutputStream(new File(outputLocation.getAbsolutePath() + File.separator + filename));
                response.setContentType("text/plain");
                String exportURL = gsver.getRequest().getContextPath() + "/" + relativeOutputFolder.replace(File.separator, "/") + filename;
                LOG.debug("On-server export file for export " + processId + ": " + exportURL);
                response.getWriter().write(exportURL);
                response.flushBuffer();
            } else {
                os = response.getOutputStream();
                response.setContentType("application/zip");
                response.setHeader("Content-disposition", "inline; filename=" + filename);
            }

            GenotypingProject project = mongoTemplate.findById(projId, GenotypingProject.class);
			Collection<GenotypingSample> samples1 = MgdbDao.getSamplesForProject(sModule, projId, selectedIndividualList1);
			Collection<GenotypingSample> samples2 = MgdbDao.getSamplesForProject(sModule, projId, selectedIndividualList2);

            Map<String, InputStream> readyToExportFiles = new HashMap<>();
            String sCitingText = appConfig.get("howToCite");
            if (sCitingText == null)
            	sCitingText = "Please cite Gigwa as follows:\nGuilhem Sempéré, Adrien Pétel, Mathieu Rouard, Julien Frouin, Yann Hueber, Fabien De Bellis, Pierre Larmande,\nGigwa v2—Extended and improved genotype investigator, GigaScience, Volume 8, Issue 5, May 2019, giz051, https://doi.org/10.1093/gigascience/giz051";
            String projDesc = project.getDescription();
            if (projDesc != null && projDesc.contains("HOW TO CITE"))
            	sCitingText += (sCitingText.length() > 0 ? "\n\n" : "") + "Please cite project data as follows:\n" + projDesc.substring(projDesc.indexOf("HOW TO CITE") + 11).replaceAll("\n\n*", "\n").trim();
            if (sCitingText.length() > 0)
            	readyToExportFiles.put("HOW_TO_CITE.txt", new ByteArrayInputStream(sCitingText.getBytes("UTF-8")));

            final OutputStream finalOS = os;
			ArrayList<GenotypingSample> samplesToExport = MgdbDao.getSamplesForProject(sModule, projId, individualsToExport);
            if (individualOrientedExportHandler != null)
            {
                if (!progress.isAborted()) {
                    Thread exportThread = new Thread() {
                    	public void run() {
                            try {
                                progress.addStep("Reading and re-organizing genotypes"); // initial step will consist in organizing genotypes by individual rather than by marker
                                progress.moveToNextStep();	// done with identifying variants
                				Map<String, File> exportFiles = individualOrientedExportHandler.createExportFiles(sModule, markerCursor.copy(), samples1, samples2, processId, gsver.getAnnotationFieldThresholds(), gsver.getAnnotationFieldThresholds2(), samplesToExport, progress);

                				for (String step : individualOrientedExportHandler.getStepList())
                                    progress.addStep(step);
                                progress.moveToNextStep();
								individualOrientedExportHandler.exportData(finalOS, sModule, exportFiles.values(), true, progress, markerCursor, null, readyToExportFiles);
					            if (!progress.isAborted()) {
					                LOG.info("doVariantExport took " + (System.currentTimeMillis() - before) / 1000d + "s to process " + count + " variants and " + CollectionUtils.union(selectedIndividualList1, selectedIndividualList2).size() + " individuals");
					                progress.markAsComplete();
					            }
							}
                            catch (Exception e) {
					            LOG.error("Error exporting data", e);
					            progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
							}
                            finally {
                                markerCursor.close();
                                try
                                {
									finalOS.close();
								}
                                catch (IOException ignored)
                                {}
                            }
                    	}
                    };
                    if (gsver.isKeepExportOnServer())
                    	exportThread.start();
                    else
                    {
                        String contentType = individualOrientedExportHandler.getExportContentType();
                        if (contentType != null && contentType.trim().length() > 0)
                        	response.setContentType(contentType);

                        exportThread.run();
                    }
                }
            }
            else if (markerOrientedExportHandler != null)
            {
                for (String step : markerOrientedExportHandler.getStepList()) {
                    progress.addStep(step);
                }
                progress.moveToNextStep();	// done with identifying variants

                String contentType = markerOrientedExportHandler.getExportContentType();
                if (contentType != null && contentType.trim().length() > 0)
                	response.setContentType(contentType);
                
                Thread exportThread = new Thread() {
                	public void run() {
                        try {
                        	markerOrientedExportHandler.exportData(finalOS, sModule, samples1, samples2, progress, markerCursor, null, gsver.getAnnotationFieldThresholds(), gsver.getAnnotationFieldThresholds2(), samplesToExport, readyToExportFiles);
                            if (!progress.isAborted()) {
                                LOG.info("doVariantExport took " + (System.currentTimeMillis() - before) / 1000d + "s to process " + count + " variants and " + CollectionUtils.union(selectedIndividualList1, selectedIndividualList2).size() + " individuals");
                                progress.markAsComplete();
                            }
						}
                        catch (Exception e) {
				            LOG.error("Error exporting data", e);
				            progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
						}
                        finally {
                            markerCursor.close();
                            try
                            {
								finalOS.close();
							}
                            catch (IOException ignored)
                            {}
                        }
                	}
                };
                if (gsver.isKeepExportOnServer())
                	exportThread.start();
                else
                	exportThread.run();
            }
			else
				throw new Exception("No export handler found for format " + gsver.getExportFormat());
        } catch (Throwable t) {
            LOG.error("Error exporting data", t);
            progress.setError("Error exporting data: " + t.getClass().getSimpleName() + (t.getMessage() != null ? " - " + t.getMessage() : ""));
            return;
        }
    }

    @Override
    public int getSequenceFilterCount(HttpServletRequest request, String sModule) throws IOException {
        int result = -1;
        File sequenceListFile = new File(request.getSession().getServletContext().getRealPath(SEQLIST_FOLDER + File.separator + request.getSession().getId() + "_" + sModule));
        if (sequenceListFile.exists() && sequenceListFile.length() > 0) {
            try (LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(sequenceListFile))) {
                lineNumberReader.skip(Long.MAX_VALUE);
                result = lineNumberReader.getLineNumber();
            }
        }
        return result;
    }

    @Override
    public ArrayList<String> getSequenceIDsBeingFilteredOn(HttpSession session, String sModule) {
        ArrayList<String> sequences = new ArrayList<>();
        File selectionFile = new File(session.getServletContext().getRealPath(SEQLIST_FOLDER) + File.separator + session.getId() + "_" + sModule);
        if (selectionFile.exists() && selectionFile.length() > 0) {
            Scanner sc = null;
            try {
                sc = new Scanner(selectionFile);
                sc.nextLine();	// skip queryKey line
                while (sc.hasNextLine()) {
                    sequences.add(sc.nextLine().trim());

                } // skip queryKey line
            } catch (FileNotFoundException ex) {
                LOG.debug("couldn't find sequence list file", ex);
            } finally {
                sc.close();
            }
        }
        return sequences.isEmpty() ? null : sequences;
    }

    @Override
    public void clearSequenceFilterFile(HttpServletRequest request, String sModule) {
        File selectionFile = new File(request.getSession().getServletContext().getRealPath(SEQLIST_FOLDER) + File.separator + request.getSession().getId() + "_" + sModule);
        if (selectionFile.exists()) {
            selectionFile.delete();
        }
    }

    @Override
    public void cleanupExpiredExportData(HttpServletRequest request) throws IOException {
        if (request.getSession() == null) {
            return;	// working around some random bug
        }
        long nowMillis = new Date().getTime();
        File filterOutputLocation = new File(request.getSession().getServletContext().getRealPath(FRONTEND_URL + File.separator + TMP_OUTPUT_FOLDER));
        if (filterOutputLocation.exists() && filterOutputLocation.isDirectory()) {
            for (File f : filterOutputLocation.listFiles()) {
                if (f.isDirectory() && nowMillis - f.lastModified() > EXPORT_EXPIRATION_DELAY_MILLIS) {
                    FileUtils.deleteDirectory(f);	// it is an expired job-output-folder
                    LOG.info("Temporary folder was deleted: " + f.getPath());
                }
            }
        }
    }

    @Override
    public boolean abortProcess(String processID) {
        ProgressIndicator progress = ProgressIndicator.get(processID);
        if (progress != null) {
            progress.abort();
            LOG.debug("Aborting process: " + processID + " [" + progress.hashCode() + "]");
            return true;
        }
        return false;
    }

    @Override
    public void onInterfaceUnload(String sModule, String processID) {
        String collName = getTemporaryVariantCollection(sModule, processID, false).getName();
        MongoTemplateManager.get(sModule).dropCollection(collName);        
        LOG.debug("Dropped collection " + sModule + "." + collName);
    }

    @Override
    public Collection<String> distinctSequencesInSelection(HttpServletRequest request, String sModule, int projId, String processID) {
        String sShortProcessID = processID/*.substring(1 + processID.indexOf('|'))*/;
        DBCollection tmpVarColl = getTemporaryVariantCollection(sModule, sShortProcessID, false);
        if (tmpVarColl.count() == 0) {
            return listSequences(request, sModule, projId);	// working on full dataset
        }
        List<String> distinctSequences = tmpVarColl.distinct(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE);
        TreeSet<String> sortedResult = new TreeSet<>(new AlphaNumericComparator());
        sortedResult.addAll(distinctSequences);
        return sortedResult;
    }

    @Override
    public String getSequenceFilterQueryKey(HttpServletRequest request, String sModule) throws IOException {

        String qk = null;
        File sequenceListFile = new File(request.getSession().getServletContext().getRealPath(SEQLIST_FOLDER + File.separator + request.getSession().getId() + "_" + sModule));
        if (sequenceListFile.exists() && sequenceListFile.length() > 0) {
            try (LineNumberReader lineNumberReader = new LineNumberReader(new FileReader(sequenceListFile))) {
                qk = lineNumberReader.readLine();
            }
        }
        return qk;
    }

    @Override
    public String getQueryKey(GigwaSearchVariantsRequest gsvr) {
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        int projId = Integer.parseInt(info[1]);
    	String queryKey = projId + ":"
    			+ gsvr.getSelectedVariantTypes() + ":"
    			+ gsvr.getReferenceName() + ":"
    			+ (gsvr.getStart() == null ? "" : gsvr.getStart()) + ":"
    			+ (gsvr.getEnd() == null ? "" : gsvr.getEnd()) + ":"
    			+ gsvr.getAlleleCount() + ":"
    			+ gsvr.getGeneName() + ":"
    			
    			+ gsvr.getCallSetIds() + ":"
    			+ gsvr.getAnnotationFieldThresholds() + ":"
    			+ gsvr.getGtPattern() + ":"
    			+ gsvr.getMostSameRatio() + ":"
    			+ gsvr.getMissingData() + ":"
    			+ gsvr.getMinmaf() + ":"
    			+ gsvr.getMaxmaf() + ":"
    			
    			+ gsvr.getCallSetIds2() + ":"
    			+ gsvr.getAnnotationFieldThresholds2() + ":"
    			+ gsvr.getGtPattern2() + ":"
    			+ gsvr.getMostSameRatio2() + ":"
    			+ gsvr.getMissingData2() + ":"
    			+ gsvr.getMinmaf2() + ":"
    			+ gsvr.getMaxmaf2() + ":"

    			+ gsvr.isDiscriminate() + ":"
    			+ gsvr.getVariantEffect();
//    	System.out.println(queryKey);
        return Helper.convertToMD5(queryKey);
	}
    
    public boolean findDefautRangeMinMax(GigwaDensityRequest gsvdr, String collectionName, ProgressIndicator progress)
	{
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvdr.getVariantSetId(), 2);
        String sModule = info[0];
        
		final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		
		BasicDBList matchAndList = new BasicDBList();
		matchAndList.add(new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, gsvdr.getDisplayedSequence()));
		if (gsvdr.getDisplayedVariantType() != null)
			matchAndList.add(new BasicDBObject(VariantData.FIELDNAME_TYPE, gsvdr.getDisplayedVariantType()));
		BasicDBObject match = new BasicDBObject("$match", new BasicDBObject("$and", matchAndList));

		BasicDBObject groupFields = new BasicDBObject("_id", null);
		groupFields.put("min", new BasicDBObject("$min", "$" + (VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE)));
		groupFields.put("max", new BasicDBObject("$max", "$" + (VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE)));
		BasicDBObject group = new BasicDBObject("$group", groupFields);

		List<DBObject> pipeline = new ArrayList<DBObject>();
		pipeline.add(match);
		pipeline.add(group);
		Iterator<DBObject> iterator = mongoTemplate.getCollection(collectionName).aggregate(pipeline).results().iterator();
		if (!iterator.hasNext())
		{
			if (progress != null)
				progress.markAsComplete();
			return false;	// no variants found matching filter
		}

		DBObject aggResult = (DBObject) iterator.next();
		if (gsvdr.getDisplayedRangeMin() == null)
			gsvdr.setDisplayedRangeMin((Long) aggResult.get("min"));
		if (gsvdr.getDisplayedRangeMax() == null)
			gsvdr.setDisplayedRangeMax((Long) aggResult.get("max"));
		return true;
	}

    @Override
    public Map<Long, Long> selectionDensity(GigwaDensityRequest gdr) throws Exception {
		long before = System.currentTimeMillis();

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gdr.getVariantSetId(), 2);
        String sModule = info[0];
        
		ProgressIndicator progress = new ProgressIndicator(tokenManager.readToken(gdr.getRequest()), new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "variant density on sequence " + gdr.getDisplayedSequence()});
		ProgressIndicator.registerProgressIndicator(progress);

		final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        final BasicDBList variantQueryDBList = (BasicDBList) buildVariantDataQuery(gdr, getSequenceIDsBeingFilteredOn(gdr.getRequest().getSession(), sModule));
		
		DBCollection tmpVarColl = getTemporaryVariantCollection(sModule, tokenManager.readToken(gdr.getRequest()), false);
		long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getName());
		if (GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr).size() > 0 && nTempVarCount == 0)
		{
			progress.setError(MESSAGE_TEMP_RECORDS_NOT_FOUND);
			return null;
		}

		final String mainVarCollName = mongoTemplate.getCollectionName(VariantData.class), usedVarCollName = nTempVarCount == 0 ? mainVarCollName : tmpVarColl.getName();
		final ConcurrentHashMap<Long, Long> result = new ConcurrentHashMap<Long, Long>();

		if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
			if (!findDefautRangeMinMax(gdr, usedVarCollName, progress))
				return result;

		final AtomicInteger nTotalTreatedVariantCount = new AtomicInteger(0);
		final int intervalSize = Math.max(1, (int) ((gdr.getDisplayedRangeMax() - gdr.getDisplayedRangeMin()) / gdr.getDisplayedRangeIntervalCount()));
		final ArrayList<Thread> threadsToWaitFor = new ArrayList<Thread>();
		final long rangeMin = gdr.getDisplayedRangeMin();
		final ProgressIndicator finalProgress = progress;

		for (int i=0; i<gdr.getDisplayedRangeIntervalCount(); i++)
		{
			BasicDBList queryList = new BasicDBList();
			queryList.add(new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, gdr.getDisplayedSequence()));
			if (gdr.getDisplayedVariantType() != null)
				queryList.add(new BasicDBObject(VariantData.FIELDNAME_TYPE, gdr.getDisplayedVariantType()));
			String startSitePath = VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE;
			queryList.add(new BasicDBObject(startSitePath, new BasicDBObject("$gte", gdr.getDisplayedRangeMin() + (i*intervalSize))));
			queryList.add(new BasicDBObject(startSitePath, new BasicDBObject(i < gdr.getDisplayedRangeIntervalCount() - 1 ? "$lt" : "$lte", i < gdr.getDisplayedRangeIntervalCount() - 1 ? gdr.getDisplayedRangeMin() + ((i+1)*intervalSize) : gdr.getDisplayedRangeMax())));
			if (nTempVarCount == 0 && !variantQueryDBList.isEmpty())
				queryList.addAll(variantQueryDBList);
			final long chunkIndex = i;

            Thread t = new Thread() {
            	public void run() {
            		if (!finalProgress.isAborted())
            		{
	        			long partialCount = mongoTemplate.getCollection(usedVarCollName).count(new BasicDBObject("$and", queryList));
	        			nTotalTreatedVariantCount.addAndGet((int) partialCount);
	        			result.put(rangeMin + (chunkIndex*intervalSize), partialCount);
	        			finalProgress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
            		}
            	}
            };

            if (chunkIndex%INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS  == (INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS-1))
            	t.run();	// run synchronously
            else
            {
            	threadsToWaitFor.add(t);
            	t.start();	// run asynchronously for better speed
            }
		}

		if (progress.isAborted())
			return null;

		for (Thread ttwf : threadsToWaitFor)	// wait for all threads before moving to next phase
			ttwf.join();

		progress.setCurrentStepProgress(100);
		LOG.debug("selectionDensity treated " + nTotalTreatedVariantCount.get() + " variants in " + (System.currentTimeMillis() - before)/1000f + "s");
		progress.markAsComplete();

		return new TreeMap<Long, Long>(result);
	}

    @Override
    public Map<Long, Integer> selectionVcfFieldPlotData(GigwaVcfFieldPlotRequest gvfpr) throws Exception {
		long before = System.currentTimeMillis();

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gvfpr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);
        
		ProgressIndicator progress = new ProgressIndicator(tokenManager.readToken(gvfpr.getRequest()), new String[] {"Calculating plot data for " + gvfpr.getVcfField() +  " field regarding " + (gvfpr.getDisplayedVariantType() != null ? gvfpr.getDisplayedVariantType() + " " : "") + "variants on sequence " + gvfpr.getDisplayedSequence()});
		ProgressIndicator.registerProgressIndicator(progress);

		final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        final BasicDBList variantQueryDBList = (BasicDBList) buildVariantDataQuery(gvfpr, getSequenceIDsBeingFilteredOn(gvfpr.getRequest().getSession(), sModule));

		DBCollection tmpVarColl = getTemporaryVariantCollection(sModule, tokenManager.readToken(gvfpr.getRequest()), false);
		long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getName());
		if (GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gvfpr).size() > 0 && nTempVarCount == 0)
		{
			progress.setError(MESSAGE_TEMP_RECORDS_NOT_FOUND);
			return null;
		}

		final String mainVarCollName = mongoTemplate.getCollectionName(VariantData.class), usedVarCollName = nTempVarCount == 0 ? mainVarCollName : tmpVarColl.getName();
		final ConcurrentHashMap<Long, Integer> result = new ConcurrentHashMap<Long, Integer>();

		if (gvfpr.getDisplayedRangeMin() == null || gvfpr.getDisplayedRangeMax() == null)
			if (!findDefautRangeMinMax(gvfpr, usedVarCollName, progress))
				return result;

		final int intervalSize = Math.max(1, (int) ((gvfpr.getDisplayedRangeMax() - gvfpr.getDisplayedRangeMin()) / gvfpr.getDisplayedRangeIntervalCount()));
		final ArrayList<Thread> threadsToWaitFor = new ArrayList<Thread>();
		final long rangeMin = gvfpr.getDisplayedRangeMin();
		final ProgressIndicator finalProgress = progress;
			
		HashMap<Integer, List<Integer>> individualIndexToSampleListMap = new HashMap<Integer, List<Integer>>();
		if (gvfpr.getPlotIndividuals() == null || gvfpr.getPlotIndividuals().size() == 0)
			gvfpr.setPlotIndividuals(MgdbDao.getProjectIndividuals(sModule, projId));
		for (int k=0; k<gvfpr.getPlotIndividuals().size(); k++)
		{
			String ind = gvfpr.getPlotIndividuals().get(k);
			List<Integer> sampleIndexes = MgdbDao.getSamplesForProject(sModule, projId, Arrays.asList(ind)).stream().map(sp -> sp.getId()).collect(Collectors.toList());
			individualIndexToSampleListMap.put(k, sampleIndexes);
		}
		
		for (int i=0; i<gvfpr.getDisplayedRangeIntervalCount(); i++)
		{
			List<Criteria> crits = new ArrayList<Criteria>();
			crits.add(Criteria.where(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE).is(gvfpr.getDisplayedSequence()));
			if (gvfpr.getDisplayedVariantType() != null)
				crits.add(Criteria.where(VariantData.FIELDNAME_TYPE).is(gvfpr.getDisplayedVariantType()));
			String startSitePath = VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE;
			crits.add(Criteria.where(startSitePath).gte(gvfpr.getDisplayedRangeMin() + (i*intervalSize)));
			if (i < gvfpr.getDisplayedRangeIntervalCount() - 1)
				crits.add(Criteria.where(startSitePath).lt(gvfpr.getDisplayedRangeMin() + ((i+1)*intervalSize)));
			else
				crits.add(Criteria.where(startSitePath).lte(gvfpr.getDisplayedRangeMax()));

			final Query query = new Query(new Criteria().andOperator(crits.toArray(new Criteria[crits.size()])));
						
			final long chunkIndex = i;
            Thread t = new Thread() {
            	public void run() {
            		if (!finalProgress.isAborted())
            		{
            			List<Comparable> variantsInInterval = mongoTemplate.getCollection(usedVarCollName).distinct("_id", query.getQueryObject());

            			final ArrayList<BasicDBObject> pipeline = new ArrayList<BasicDBObject>();

            			BasicDBObject group = new BasicDBObject();
            			ArrayList<Object> individualValuePaths = new ArrayList<>();
            			for (int j=0; j<gvfpr.getPlotIndividuals().size(); j++)
            			{
            				List<Integer> individualSamples = individualIndexToSampleListMap.get(j);
            				if (individualSamples.size() == 1)
            					individualValuePaths.add("$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + individualSamples.get(0) + "." + SampleGenotype.SECTION_ADDITIONAL_INFO + "." + gvfpr.getVcfField());
            				else
            				{	// only take into account the sample with the highest value
            					if (group.size() == 0)
            						group.put("_id", null);
            					for (int l=0; l<variantsInInterval.size(); l++)
            					{
                					ArrayList<BasicDBObject> sampleFields = new ArrayList<>();
	            					for (int k=0; k<individualSamples.size(); k++)
	            					{
	            						if (l == 0)
	            							group.put(gvfpr.getVcfField() + individualSamples.get(k), new BasicDBObject("$addToSet", "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + individualSamples.get(k) + "." + SampleGenotype.SECTION_ADDITIONAL_INFO + "." + gvfpr.getVcfField()));
	                					sampleFields.add(new BasicDBObject("$arrayElemAt", new Object[] {"$" + gvfpr.getVcfField() + individualSamples.get(k), l}));
	            					}
                					individualValuePaths.add(new BasicDBObject("$max", sampleFields));
            					}
            				}
            			}
            			if (group.size() > 0)
            				pipeline.add(new BasicDBObject("$group", group));
            			pipeline.add(new BasicDBObject("$project", new BasicDBObject(gvfpr.getVcfField(), new BasicDBObject("$sum", individualValuePaths))));
            			
            			BasicDBList matchList = new BasicDBList();
            			matchList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", variantsInInterval)));
            			if (nTempVarCount == 0 && !variantQueryDBList.isEmpty())
            				matchList.addAll(variantQueryDBList);
            			pipeline.add(0, new BasicDBObject("$match", new BasicDBObject("$and", matchList)));
            			
	        			Iterator<DBObject> it = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).aggregate(pipeline).results().iterator();
	        			result.put(rangeMin + (chunkIndex*intervalSize), it.hasNext() ? Double.valueOf(it.next().get(gvfpr.getVcfField()).toString()).intValue() : 0);
	        			finalProgress.setCurrentStepProgress((short) result.size() * 100 / gvfpr.getDisplayedRangeIntervalCount());
            		}
            	}
            };

            if (chunkIndex%INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS  == (INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS-1))
            	t.run();	// run synchronously
            else
            {
            	threadsToWaitFor.add(t);
            	t.start();	// run asynchronously for better speed
            }
		}

		if (progress.isAborted())
			return null;

		for (Thread ttwf : threadsToWaitFor)	// wait for all threads before moving to next phase
			ttwf.join();

		progress.setCurrentStepProgress(100);
		LOG.debug("selectionVcfFieldPlotData treated " + gvfpr.getDisplayedRangeIntervalCount() + " intervals on sequence " + gvfpr.getDisplayedSequence() + " between " + gvfpr.getDisplayedRangeMin() + " and " + gvfpr.getDisplayedRangeMax() + " in " + (System.currentTimeMillis() - before)/1000f + "s");
		progress.markAsComplete();

		return new TreeMap<Long, Integer>(result);
	}

    @Override
    public List<String> getProjectSequences(String sModule, int projId) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        GenotypingProject proj = mongoTemplate.findById(projId, GenotypingProject.class);
        return new ArrayList<String>(proj.getSequences());
    }

    /**
     * Import a sequence from a fasta in the database
     *
     * @param module
     * @param filePath fasta file local or remote path (accept ftp:// )
     * @param mode
     * @return boolean true if import succeded
     * @throws java.lang.Exception
     */
    public boolean importSequenceInDB(String module, String filePath, String mode) throws Exception {
        boolean success = false;
        try {
            SequenceImport.main(new String[]{module, filePath, mode});
            success = true;
        } catch (Exception ex) {
            throw ex;
        }
        return success;
    }

    /**
     * get all run in a project
     *
     * @param module
     * @param projId
     * @return
     */
    public List<String> getRunList(String module, int projId) {

        List<String> listRun;
        MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_RUNS);
       	q.addCriteria(Criteria.where("_id").is(projId));
        GenotypingProject project = mongoTemplate.findOne(q, GenotypingProject.class);
        listRun = project.getRuns();
        if (listRun == null) {
            return new ArrayList<>();
        } else {
            return listRun;
        }
    }

    /**
     * get information on a list of sequence
     *
     * @param module
     * @param sequenceList
     * @return
     */
    public Map<String, Map<String, Object>> getSequenceInfo(String module, List<String> sequenceList) {

        Map<String, Map<String, Object>> listSeqInfo = new HashMap<>();
        Map<String, Object> info = new HashMap<>();
        DBObject sequence;

        // no need to filter on projId since we are searching on sequenceId 
        ArrayList<DBObject> aggregationParam = new ArrayList<>();
        aggregationParam.add(new BasicDBObject("$match", new BasicDBObject("_id", new BasicDBObject("$in", sequenceList))));

        Cursor genotypingDataCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(aggregationParam, AggregationOptions.builder().allowDiskUse(isAggregationAllowedToUseDisk()).build());

        if (genotypingDataCursor != null && genotypingDataCursor.hasNext()) {

            while (genotypingDataCursor.hasNext()) {

                sequence = genotypingDataCursor.next();

                // empty the table
                info.clear();

                // no need for sequence? 
                info.put("sequence", sequence.get(Sequence.FIELDNAME_SEQUENCE));

                info.put("length", sequence.get(Sequence.FIELDNAME_LENGTH));
                info.put("cheksum", sequence.get(Sequence.FIELDNAME_CHECKSUM));

                listSeqInfo.put((String) sequence.get("_id"), info);

            }
        }
        return listSeqInfo;
    }

	/**
	 * get a list of variant in ga4gh format from a DBCursor
	 *
	 * @param module
	 * @param projId
	 * @param cursor
	 * @param samples
	 * @param run
	 * @return List<Variant>
	 * @throws AvroRemoteException
	 */
	private List<Variant> getVariantListFromDBCursor(String module, int projId, DBCursor cursor, Collection<GenotypingSample> samples, String run) throws AvroRemoteException
	{
//    	long before = System.currentTimeMillis();
        LinkedHashMap<Comparable, Variant> varMap = new LinkedHashMap<>();
        LinkedHashMap<Comparable, String> typeMap = new LinkedHashMap<>();

        // parse the cursor to create all GAVariant 
        while (cursor.hasNext()) {

            DBObject obj = cursor.next();
            // save the Id of each variant in the cursor 
            String id = (String) obj.get("_id");
            List<String> knownAlleles = ((List<String>) obj.get(VariantData.FIELDNAME_KNOWN_ALLELE_LIST));

            DBObject rp = ((DBObject) obj.get(VariantData.FIELDNAME_REFERENCE_POSITION));
            if (rp == null)
            	throw new AvroRemoteException("Variant " + id.toString() + " has no reference position!");
            
            Long start = (Long) rp.get(ReferencePosition.FIELDNAME_START_SITE);
            Long end = (Long) rp.get(ReferencePosition.FIELDNAME_END_SITE);
            if (end == null && start != null)
            	end = start;

            typeMap.put(id, (String) obj.get(VariantData.FIELDNAME_TYPE));
            Variant.Builder variantBuilder = Variant.newBuilder()
                    .setId(createId(module, projId, id.toString()))
                    .setVariantSetId(createId(module, projId))
                    .setReferenceName((String) rp.get(ReferencePosition.FIELDNAME_SEQUENCE))
                    .setStart(start)
                    .setEnd(end);
            if (knownAlleles.size() == 0)
            	throw new AvroRemoteException("Variant " + id.toString() + " has no known alleles!");

            variantBuilder.setReferenceBases(knownAlleles.get(0)); // reference is the first one in VCF files
           	variantBuilder.setAlternateBases(knownAlleles.subList(1, knownAlleles.size()));
           	varMap.put(id, variantBuilder.build());
        }

        // get the VariantRunData containing annotations
        ArrayList<DBObject> pipeline = new ArrayList<>();
        // wanted fields 
        BasicDBObject fields = new BasicDBObject();
        fields.put(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE, 1);
        fields.put(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, 1);

        // get the genotype for wanted individuals/callSet only 
		final Map<Integer, String> sampleIdToIndividualMap = new HashMap<>();
        for (GenotypingSample sample : samples)
	        if (run == null || sample.getRun().equals(run)) {
	            fields.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + sample.getId(), 1);
	            sampleIdToIndividualMap.put(sample.getId(), sample.getIndividual());
	        }
        
        BasicDBList matchAndList = new BasicDBList();
        matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", varMap.keySet())));
        matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID, projId));
        if (run != null)
        	matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_RUNNAME, run));
        pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", matchAndList)));
        pipeline.add(new BasicDBObject("$project", fields));
        
        HashSet<String> variantsForWhichAnnotationWasRetrieved = new HashSet<>();

        Cursor genotypingDataCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(pipeline, AggregationOptions.builder().allowDiskUse(isAggregationAllowedToUseDisk()).build());
        while (genotypingDataCursor.hasNext())
        {
            DBObject variantObj = genotypingDataCursor.next();
            String varId = (String) Helper.readPossiblyNestedField(variantObj, "_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
        	Variant var = varMap.get(varId);
        	if (var == null /* should not happen! */|| variantsForWhichAnnotationWasRetrieved.contains(varId))
        		continue;
        	
        	variantsForWhichAnnotationWasRetrieved.add(varId);
            TreeSet<Call> calls = new TreeSet(new AlphaNumericComparator<Call>());	// for automatic sorting

            Map<String, List<String>> annotations = new HashMap<>();
            List<String> infoType = new ArrayList<>();
            infoType.add(typeMap.get(varId));
            annotations.put("type", infoType);

            // for each annotation field
            for (String key : variantObj.keySet()) {

                switch (key) {

                    // this goes in Call  || should not be called if sp field is not present 
                    case VariantRunData.FIELDNAME_SAMPLEGENOTYPES:
                        // get genotype map 
                        Map<String, Object> callMap = (Map<String, Object>) variantObj.get(key);

                        // for each individual/CallSet
                        for (Integer sampleId : sampleIdToIndividualMap.keySet()) {
                            DBObject callObj = (DBObject) callMap.get("" + sampleId);
                            double[] gl;
                            List<Double> listGL = new ArrayList<>();

                            List<Integer> genotype = new ArrayList<>();
                            Map<String, List<String>> aiCall = new HashMap<>();
                            String phaseSet = null;
                            
                            if (callObj != null)
                            {
                                Map<String, Object> callAdditionalInfo = (Map<String, Object>) callObj.get("ai");

                                // if field ai is present 
                                if (callAdditionalInfo != null)
                                    for (String aiKey : callAdditionalInfo.keySet()) {
                                    	if (aiKey.equals(VCFConstants.GENOTYPE_PL_KEY))
                                    	{
                                            gl = GenotypeLikelihoods.fromPLField(callAdditionalInfo.get(aiKey).toString()).getAsVector();
                                            for (int h = 0; h < gl.length; h++)
                                                listGL.add(gl[h]);
                                    	}
                                        switch (aiKey) {
                                            case VCFConstants.GENOTYPE_LIKELIHOODS_KEY:
                                                listGL = (List<Double>) callAdditionalInfo.get(aiKey);
                                                break;

                                            case VCFConstants.PHASE_SET_KEY:
                                            case VariantData.GT_FIELD_PHASED_ID:
                                                phaseSet = callAdditionalInfo.get(aiKey).toString();
                                                break;

                                            default:
                                                List<String> callAddinfoContent = new ArrayList<>();
                                                callAddinfoContent.add(callAdditionalInfo.get(aiKey).toString());
                                                aiCall.put(aiKey, callAddinfoContent);
                                                break;
                                        }
                                    }

                                // get GT info 
                                String gt = (String) callObj.get("gt");

                                if (gt == null || gt.startsWith(".")) {
                                    // if we don't know the genotype, do nothing
                                } else {
                                    String[] gen;
                                    if (gt.contains("/")) {
                                        gen = gt.split("/");
                                    } else {
                                        gen = gt.split(GigwaGa4ghServiceImpl.ID_SEPARATOR);
                                    }
                                    for (String gen1 : gen) {
                                        genotype.add(Integer.parseInt(gen1));
                                    }
                                }
                            }
                            Call call = Call.newBuilder()
                                    .setCallSetId(createId(module, projId, sampleIdToIndividualMap.get(sampleId)))
                                    .setGenotype(genotype)
                                    .setGenotypeLikelihood(listGL)
                                    .setPhaseset(phaseSet)
                                    .setInfo(aiCall)
                                    .build();

                            calls.add(call);
                        }
                        break;

                    case VariantRunData.SECTION_ADDITIONAL_INFO:

                        Map<String, Object> additionalInfos = (Map<String, Object>) variantObj.get(key);

                        for (String subKey : additionalInfos.keySet()) {

                            if (subKey.equals("") || subKey.equals(VcfImport.ANNOTATION_FIELDNAME_ANN) || subKey.equals(VcfImport.ANNOTATION_FIELDNAME_CSQ)) {
                                // if VCF has empty field (";") do not retrieve it 

                                // field EFF should be stored in variantAnnotation !
                                // stored in ai for the moment, not supported by ga4gh  
                                // ANN (vcf 4.2) is stored in variantAnnotation 
                            } else if (subKey.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE)) {
                                List<String> listGene = (List<String>) additionalInfos.get(subKey);
                                annotations.put(subKey, listGene);
                            } else if (subKey.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME)) {
                                List<String> listEffect = (List<String>) additionalInfos.get(subKey);
                                annotations.put(subKey, listEffect);
                            } else {

                            }
                        }
                        break;
                    default:
                        // "_id" and "_class", do nothing 
                        break;
                }
            }
            // add the call list
            var.setCalls(new ArrayList<Call>(calls));
            // add the annotation map to the variant 
            var.setInfo(annotations);
        }

//        LOG.debug("getVariantListFromDBCursor took " + (System.currentTimeMillis() - before) / 1000f + "s for " + listVar.size() + " variants and " + samples.size() + " samples");
        return new ArrayList<Variant>(varMap.values());
    }

    /**
     * create ID from a list of param
     *
     * @param param
     * @return String the id
     */
    public String createId(Comparable... params) {
        String result = "";

        for (Comparable val : params) {
            result += val + ID_SEPARATOR;
        }
        return result.substring(0, result.length() - 1);
    }

    /**
     * get the Metadata number from a VCFCompoundLine
     *
     * @param vcf
     * @return String number
     */
    private String getNumber(VCFCompoundHeaderLine vcf) {

        String number;
        switch (vcf.getCountType()) {
            case A:
                number = "";
                break;
            case G:
                number = "";
                break;
            case INTEGER:
                number = Integer.toString(vcf.getCount());
                break;
            case R:
                number = "";
                break;
            case UNBOUNDED:

                // ga4gh python serveur return "." when unbounded 
                // but no information about it in ga4gh documentation 
                number = Integer.toString(-1);
                break;
            default:
                number = "";
                break;
        }
        return number;
    }

    /**
     * return the progress
     *
     * @param processID
     * @return ProgressIndicator
     */
    public ProgressIndicator getProgressIndicator(String processID) {
        return ProgressIndicator.get(processID/*.substring(1 + processID.indexOf('|'))*/);
    }

    /**
     * get the metadata for a Project/VariantSet
     *
     * @param module
     * @param listProj
     * @return List<List<VariantSetMetadata>>
     */
    private List<VariantSetMetadata> getMetadataList(String module, String proj) {

        List<VariantSetMetadata> listMetadata;
        VariantSetMetadata metadata;
        Map<String, VariantSetMetadata> mapMetadata = new HashMap<>();
        DBVCFHeader vcfHeader;

        VCFInfoHeaderLine vci;
        VCFFilterHeaderLine vcf;
        VCFFormatHeaderLine vcfo;
        VCFSimpleHeaderLine vcm;
        VCFHeaderLine vcmo;
        String metadataKey;
        int i;
        Map<String, List<String>> info = new HashMap<>();

        // get the list of vcf header for this project         
        BasicDBObject whereQuery = new BasicDBObject();
        whereQuery.put("_id." + DBVCFHeader.VcfHeaderId.FIELDNAME_PROJECT, Integer.parseInt(proj));
        DBCursor cursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).find(whereQuery);

        if (cursor != null && cursor.hasNext()) {

            while (cursor.hasNext()) {

                // get the vcf header of each run of a project 
                vcfHeader = DBVCFHeader.fromDBObject(cursor.next());

                // used to create id for each row header
                i = 0;

                // fill metadata list 
                for (String key : vcfHeader.getmInfoMetaData().keySet()) {

                    vci = vcfHeader.getmInfoMetaData().get(key);

                    // store the key to make sure a header is only present once 
                    // (if project has multiple run) 
                    metadataKey = vci.getKey() + "." + vci.getID();

                    metadata = VariantSetMetadata.newBuilder()
                            .setDescription(vci.getDescription())
                            .setType(vci.getType().toString())
                            .setValue(vci.getValue())
                            .setId(Integer.toString(i))
                            .setKey(metadataKey)
                            .setNumber(getNumber(vci))
                            .setInfo(info)
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);
                }
                for (String key : vcfHeader.getmFilterMetaData().keySet()) {

                    vcf = vcfHeader.getmFilterMetaData().get(key);
                    metadataKey = vcf.getKey() + "." + vcf.getID();

                    metadata = VariantSetMetadata.newBuilder()
                            .setId(Integer.toString(i))
                            .setType("String")
                            .setKey(metadataKey)
                            .setValue(vcf.getValue())
                            .setInfo(info)
                            .setDescription("")
                            .setNumber("")
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);

                }
                for (String key : vcfHeader.getmFormatMetaData().keySet()) {

                    vcfo = vcfHeader.getmFormatMetaData().get(key);

                    metadataKey = vcfo.getKey() + "." + vcfo.getID();

                    metadata = VariantSetMetadata.newBuilder()
                            .setId(Integer.toString(i))
                            .setKey(metadataKey)
                            .setValue(vcfo.getValue())
                            .setType(vcfo.getType().toString())
                            .setDescription(vcfo.getDescription())
                            .setInfo(info)
                            .setNumber(getNumber(vcfo))
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);

                }
                for (String key : vcfHeader.getmMetaData().keySet()) {

                    vcm = vcfHeader.getmMetaData().get(key);
                    metadataKey = vcm.getKey() + "." + vcm.getID();

                    metadata = VariantSetMetadata.newBuilder()
                            .setId(Integer.toString(i))
                            .setKey(metadataKey)
                            .setValue(vcm.getValue())
                            .setDescription("")
                            .setNumber("")
                            .setType("String")
                            .setInfo(info)
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);

                }

                for (String key : vcfHeader.getmOtherMetaData().keySet()) {

                    vcmo = vcfHeader.getmOtherMetaData().get(key);
                    metadataKey = vcmo.getKey();

                    metadata = VariantSetMetadata.newBuilder()
                            .setId(Integer.toString(i))
                            .setKey(metadataKey)
                            .setValue(vcmo.getValue())
                            .setDescription("")
                            .setType("String")
                            .setNumber("")
                            .setInfo(info)
                            .build();
                    i++;
                    mapMetadata.put(metadataKey, metadata);
                }
            }
            cursor.close();
        }
        // get a list of variantSetMetadata from the map
        listMetadata = new ArrayList<>(mapMetadata.values());

        return listMetadata;
    }

    /*
     * GA4GH methods, from methods interface v0.6.1 - opencb/ga4gh 04/2016 
     */
    @Override
    public VariantSet getVariantSet(String id) throws AvroRemoteException, GAException
    {
        VariantSet variantSet = null;
        // get information from id
        String[] info = GigwaSearchVariantsRequest.getInfoFromId(id, 2);
        if (info != null)
        	try
	        {
	            String module = info[0];
	            String projId = info[1];
	            
	            Query q = new Query(Criteria.where("_id").is(Integer.parseInt(projId)));
	            q.fields().include(GenotypingProject.FIELDNAME_NAME);

	            GenotypingProject proj = MongoTemplateManager.get(module).findOne(q, GenotypingProject.class);
            	List<VariantSetMetadata> metadata = getMetadataList(module, projId);
                if (proj.getDescription() != null)
                {
                	VariantSetMetadata vsmd = new VariantSetMetadata();
                	vsmd.setKey("description");
                	vsmd.setValue(proj.getDescription());
                	metadata.add(vsmd);
                }
                // Checksum : MD5 of the upper-case sequence excluding all whitespace characters 
                // length == 0 since we don't have this information in VCF files
                variantSet = VariantSet.newBuilder()
		                .setId(id)
		                .setDatasetId(module)
		                .setName(proj.getName())
		                .setReferenceSetId(module)
		                .setMetadata(metadata)
		                .build();
	        }
        	catch (NumberFormatException nfe)
        	{}
        return variantSet;
    }

    @Override
    public Variant getVariant(String id) throws AvroRemoteException, GAException {
        return getVariantWithGenotypes(id, new ArrayList() /* all individuals */);
    }

    public Variant getVariantWithGenotypes(String id, Collection<String> listInd) throws NumberFormatException, AvroRemoteException {
        Variant variant = null;
        // get information from id
        String[] info = id.split(GigwaGa4ghServiceImpl.ID_SEPARATOR);
        if (info.length <= 2) {
            // wrong number of param or wrong module name 
        } else {
            String module = info[0];
            int projId = Integer.parseInt(info[1]);
            String name = info[2];
            String run = null;
            if (info.length == 4)
                run = info[3];

            MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
            
            // get the collection of variant 
            DBCollection variantColl = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantData.class));

//            // we need to get callSet name and position in the callSet list to get corresponding genotype 
//            GenotypingProject project = mongoTemplate.findById(projId, GenotypingProject.class);

            BasicDBObject whereQuery = new BasicDBObject();
            whereQuery.put("_id", name);
            DBCursor cursor = variantColl.find(whereQuery);

            // if there is no result, return null
            if (cursor != null && cursor.hasNext()) {
				variant = getVariantListFromDBCursor(module, Integer.parseInt(info[1]), cursor, MgdbDao.getSamplesForProject(module, projId, listInd), run).get(0);
				cursor.close();
            }
        }
        return variant;
    }

    @Override
    public CallSet getCallSet(String id) throws AvroRemoteException, GAException {

        CallSet callSet = null;

        // get information from id
        String[] info = GigwaSearchVariantsRequest.getInfoFromId(id, 3);
        if (info == null) {

            // wrong number of param or wrong module name 
        } else {
            String module = info[0];
            int projId = Integer.parseInt(info[1]);
            String name = info[2];

            List<String> listVariantSetId = new ArrayList<>();
            listVariantSetId.add(createId(module, info[1]));

            // get the callSet list of a project 
            List<String> listCallSet = MgdbDao.getProjectIndividuals(module, projId);

            // check if the callSet is in the list 
            if (listCallSet.contains(name)) {
                callSet = CallSet.newBuilder()
                        .setId(id)
                        .setName(name)
                        .setVariantSetIds(listVariantSetId)
                        .setSampleId(null)
                        .build();
            }
        }
        return callSet;
    }

    @Override
    public ReferenceSet getReferenceSet(String id) throws AvroRemoteException, GAException {

        ReferenceSet referenceSet = null;

        MongoTemplate mongoTemplate = MongoTemplateManager.get(id);
        if (mongoTemplate == null) {

        } else {
            List<String> list = new ArrayList<>();
            // get the all references of the reference Set/module 
            Cursor genotypingDataCursor = MongoTemplateManager.get(id).getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).find();
            DBObject seq;

            String concatId = "";

            if (genotypingDataCursor != null) {

                List<String> listChecksum = new ArrayList<>();

                while (genotypingDataCursor.hasNext()) {
                    seq = genotypingDataCursor.next();
                    listChecksum.add((String) seq.get(Sequence.FIELDNAME_CHECKSUM));
                }
                // sort in lexicographic order
                Collections.sort(listChecksum);
                for (String checksum : listChecksum) {
                    concatId = concatId + checksum;
                }
                genotypingDataCursor.close();
            }

            String taxon = MongoTemplateManager.getTaxonName(id);
            String species = MongoTemplateManager.getSpecies(id);
            String taxoDesc = (species != null ? "Species: " + species : "") + (taxon != null && !taxon.equals(species) ? "" : (species != null ? " ; " : "") + "Taxon: " + taxon);
            referenceSet = ReferenceSet.newBuilder()
                    .setId(id)
                    .setName(id)
                    .setMd5checksum(Helper.convertToMD5(concatId))
                    .setSourceAccessions(list)
                    .setNcbiTaxonId(MongoTemplateManager.getTaxonId(id))
                    .setDescription((taxoDesc.isEmpty() ? "" : (taxoDesc + " ; ")) + mongoTemplate.getCollection(mongoTemplate.getCollectionName(GenotypingProject.class)).distinct(GenotypingProject.FIELDNAME_SEQUENCES).size() + " references ; " + mongoTemplate.count(null, VariantData.class) + " markers")
                    .build();
        }
        return referenceSet;
    }

    @Override
    public Reference getReference(String id) throws AvroRemoteException, GAException {

        Reference reference = null;

        // get information from id
        String[] info = GigwaSearchVariantsRequest.getInfoFromId(id, 3);
        if (info == null) {

            // wrong number of param or wrong module name 
        } else {
            String module = info[0];
            int projId = Integer.parseInt(info[1]);
            String name = info[2];

            MongoTemplate mongoTemplate = MongoTemplateManager.get(module);

            GenotypingProject proj = mongoTemplate.findById(projId, GenotypingProject.class);
            Set<String> listRef = proj.getSequences();

            // check if the sequence is in the list 
            if (listRef.contains(name)) {

                long length = 0L;
                String checksum = Helper.convertToMD5("");

                Sequence sequence = mongoTemplate.findById(name, Sequence.class);
                if (sequence != null) {
                    length = sequence.getLength();
                    checksum = sequence.getChecksum();
                }

                List<String> liste = new ArrayList<>();
                // Checksum : MD5 of the upper-case sequence excluding all whitespace characters 
                // length == 0 since we don't have this information in VCF files
                reference = Reference.newBuilder().setId(id)
                        .setMd5checksum(checksum)
                        .setName(name)
                        .setLength(length)
                        .setSourceAccessions(liste)
                        .build();
            }
        }
        return reference;
    }

    @Override
    public ListReferenceBasesResponse getReferenceBases(String id, ListReferenceBasesRequest lrbr) throws AvroRemoteException, GAException
    {
           String[] info = GigwaSearchVariantsRequest.getInfoFromId(id, 3);
           String module = info[0];
           String seqName = info[2];
    	   String sequenceBases = "";
		   int spaceCode = (int) '\n';
		   int chevCode = (int) '>';
		
		   if (lrbr.getEnd() > lrbr.getStart()) {
		       Sequence seq = MongoTemplateManager.get(module).findById(seqName, Sequence.class);
		       BufferedReader br = null;
		       StringBuilder builder;
		       if ((seq != null)) {
		           try {
		               br = new BufferedReader(new InputStreamReader(new FileInputStream(seq.getFilePath())));
		               builder = new StringBuilder();
		               String line;
		               int c;
		               int pos = 0;
		               // skip line to reach sequence 
		           while ((line = br.readLine()) != null && !line.startsWith(">" + id)) {
		
		           }
		           // skip char to reach start pos
		           while ((c = br.read()) != -1 && pos < lrbr.getStart()) {
		               pos++;
		           }
		           builder.append((char) c);
		
		           while ((c = br.read()) != -1 && pos < lrbr.getEnd() && c != chevCode) {
		               if (c != spaceCode) {
		                   builder.append((char) c);
		                   pos++;
		               }
		           }
		           sequenceBases = builder.toString();
		
		       }
		       catch (IOException ex)
		       {
		           LOG.warn("could not open file : " + ex);
		       }
		       finally 
		       {
		           try 
		           {
		               if (br != null)
		                   br.close();
		           }
		           catch (IOException ex)
		           {
		               LOG.warn("could not close writer : " + ex);
		           }
		       }
		   }
		   else
		       throw new GAException("No fasta for sequence " + seqName);
		   }
		   ListReferenceBasesResponse result = new ListReferenceBasesResponse();
		   result.setSequence(sequenceBases);
		   result.setOffset(lrbr.getStart());
		   return result;
    }
    
    /* needed to be able to pass totalCount to the BrAPI v2 call */
    public class SearchCallSetsResponseWrapper extends SearchCallSetsResponse {
    	private SearchCallSetsResponse scsr;    	
    	private int totalCount;

    	public SearchCallSetsResponseWrapper(SearchCallSetsResponse scsr) {
    		this.scsr = scsr;
    	}

		public SearchCallSetsResponse getResponse() {
			return scsr;
		}

		public int getTotalCount() {
			return totalCount;
		}

		public void setTotalCount(int totalCount) {
			this.totalCount = totalCount;
		}
    }

	@Override
	public SearchCallSetsResponse searchCallSets(SearchCallSetsRequest scsr) throws AvroRemoteException, GAException {
		SearchCallSetsResponse response = null;

		Authentication auth = SecurityContextHolder.getContext().getAuthentication();
    	String sCurrentUser = auth == null || "anonymousUser".equals(auth.getName()) ? "anonymousUser" : auth.getName();
    	
		// get information from id
		String[] info = GigwaSearchVariantsRequest.getInfoFromId(scsr.getVariantSetId(), 2);
		if (info == null)
			return null;

		CallSet callSet;
		int start;
		int end;
		int pageSize;
		int pageToken = 0;
		String nextPageToken;

		String module = info[0];
		MongoTemplate mongoTemplate = MongoTemplateManager.get(module);

		List<String> listVariantSetId = new ArrayList<>();
		listVariantSetId.add(scsr.getVariantSetId());

		// build the list of individuals
		Query q = new Query(Criteria.where(GenotypingSample.FIELDNAME_PROJECT_ID).is(Integer.parseInt(info[1])));
		q.fields().include(GenotypingSample.FIELDNAME_INDIVIDUAL);
		
		Map<String, Integer> indIdToSampleIdMap = new HashMap<>();
		for (GenotypingSample sample : mongoTemplate.find(q, GenotypingSample.class))
			indIdToSampleIdMap.put(sample.getIndividual(), sample.getId());

		q = new Query(Criteria.where("_id").in(indIdToSampleIdMap.keySet()));
		q.with(new Sort(Sort.Direction.ASC, "_id"));
		long totalCount = mongoTemplate.count(q, Individual.class);
		List<Individual> listInd = mongoTemplate.find(q, Individual.class);
		q = new Query(Criteria.where("_id." + CustomIndividualMetadataId.FIELDNAME_USER).is(sCurrentUser));
		List<CustomIndividualMetadata> cimdList = mongoTemplate.find(q, CustomIndividualMetadata.class);
		if(!cimdList.isEmpty()) {
			HashMap<String /* indivID */, HashMap<String, Comparable> /* additional info */> indMetadataByIdMap = new HashMap<>();
			for (CustomIndividualMetadata cimd : cimdList)
				indMetadataByIdMap.put(cimd.getId().getIndividualId(), cimd.getAdditionalInfo());
			
			for( int i=0 ; i<listInd.size(); i++) {
				String indId = listInd.get(i).getId();
				HashMap<String, Comparable>  ai = indMetadataByIdMap.get(indId);
                if(ai != null && !ai.isEmpty())
                	listInd.get(i).getAdditionalInfo().putAll(ai);
			}
		}

		List<CallSet> listCallSet = new ArrayList<>();

		int size = listInd.size();
		// if no pageSize specified, return all results
		if (scsr.getPageSize() != null) {
			pageSize = scsr.getPageSize();
		} else {
			pageSize = size;
		}
		if (scsr.getPageToken() != null) {
			pageToken = Integer.parseInt(scsr.getPageToken());
		}

		start = pageSize * pageToken;
		if (size - start <= pageSize) {
			end = size;
			nextPageToken = null;
		} else {
			end = pageSize * (pageToken + 1);
			nextPageToken = Integer.toString(pageToken + 1);
		}

		// create a callSet for each item in the list
		for (int i = start; i < end; i++) {
			final Individual ind = listInd.get(i);
			CallSet.Builder csb = CallSet.newBuilder().setId(createId(module, info[1], ind.getId())).setName(ind.getId()).setVariantSetIds(listVariantSetId).setSampleId(createId(module, info[1], ind.getId(), indIdToSampleIdMap.get(ind.getId())));
			if (!ind.getAdditionalInfo().isEmpty())
				csb.setInfo(ind.getAdditionalInfo().keySet().stream().collect(Collectors.toMap(k -> k, k -> (List<String>) Arrays.asList(ind.getAdditionalInfo().get(k).toString()))));
			callSet = csb.build();
			listCallSet.add(callSet);
		}
		response = SearchCallSetsResponse.newBuilder().setCallSets(listCallSet).setNextPageToken(nextPageToken).build();
		
		if (!Thread.currentThread().getStackTrace()[2].getClassName().toLowerCase().contains("brapi"))
			return response;
	
		SearchCallSetsResponseWrapper wrapper = new SearchCallSetsResponseWrapper(response);
		wrapper.setTotalCount((int) totalCount);
		wrapper.setCallSets(response.getCallSets());	// so this remains compatible with ga4gh
		wrapper.setNextPageToken(response.getNextPageToken());	// so this remains compatible with ga4gh
		return wrapper;
	}
	
    @Override
    public SearchReferenceSetsResponse searchReferenceSets(SearchReferenceSetsRequest srsr) throws AvroRemoteException, GAException {
        List<String> list = new ArrayList<>();
        List<ReferenceSet> listRef = new ArrayList<>();

        List<String> listModules = new ArrayList<>(MongoTemplateManager.getAvailableModules());
        Collections.sort(listModules);
        int start;
        int end;
        int pageSize;
        int pageToken = 0;
        String nextPageToken;
        String module;

        int size = listModules.size();
        // if page size is not specified, return all results
        if (srsr.getPageSize() != null) {
            pageSize = srsr.getPageSize();
        } else {
            pageSize = size;
        }
        if (srsr.getPageToken() != null) {
            pageToken = Integer.parseInt(srsr.getPageToken());
        }

        start = pageSize * pageToken;

        // nextPageToken = null if no more result
        if (size - start <= pageSize) {
            end = size;
            nextPageToken = null;
        } else {
        	end = pageSize * (pageToken + 1);
            nextPageToken = Integer.toString(pageToken + 1);
        }

        // add a Reference Set for each existing module 
        for (int i = start; i < end; i++) {
            module = listModules.get(i);
            MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
            String taxon = MongoTemplateManager.getTaxonName(module);
            String species = MongoTemplateManager.getSpecies(module);
            String taxoDesc = (species != null ? "Species: " + species : "") + (taxon != null && !taxon.equals(species) ? (species != null ? " ; " : "") + "Taxon: " + taxon : "");
            ReferenceSet referenceSet = ReferenceSet.newBuilder()
                    .setId(module)
                    .setName(module)
                    .setMd5checksum("")	/* not supported at the time */
                    .setSourceAccessions(list)
                    .setNcbiTaxonId(MongoTemplateManager.getTaxonId(module))
                    .setDescription((taxoDesc.isEmpty() ? "" : (taxoDesc + " ; ")) + mongoTemplate.getCollection(mongoTemplate.getCollectionName(GenotypingProject.class)).distinct(GenotypingProject.FIELDNAME_SEQUENCES).size() + " references ; " + mongoTemplate.count(null, VariantData.class) + " markers")
                    .build();
            listRef.add(referenceSet);
        }

        return SearchReferenceSetsResponse.newBuilder().setReferenceSets(listRef).setNextPageToken(nextPageToken).build();
    }

    @Override
    public SearchVariantSetsResponse searchVariantSets(SearchVariantSetsRequest svsr) throws AvroRemoteException, GAException {

        SearchVariantSetsResponse response = null;

        String[] info = GigwaSearchVariantsRequest.getInfoFromId(svsr.getDatasetId(), 1);
        if (info == null)
        	return null;

        String module = info[0];
        int start;
        int end;
        int pageSize;
        int pageToken = 0;

        String nextPageToken;

        Query q = new Query();
        q.fields().include(GenotypingProject.FIELDNAME_NAME);
        q.fields().include(GenotypingProject.FIELDNAME_DESCRIPTION);
        List<GenotypingProject> listProj = MongoTemplateManager.get(module).find(q, GenotypingProject.class);
        List<VariantSet> listVariantSet = new ArrayList<>();

        int size = listProj.size();
        // if page size is not specified, return all results
        if (svsr.getPageSize() != null) {
            pageSize = svsr.getPageSize();
        } else {
            pageSize = size;
        }
        if (svsr.getPageToken() != null) {
            pageToken = Integer.parseInt(svsr.getPageToken());
        }

        start = pageSize * pageToken;
        if (size - start <= pageSize) {
            end = size;
            nextPageToken = null;
        } else {
        	end = pageSize * (pageToken + 1);
            nextPageToken = Integer.toString(pageToken + 1);
        }

        for (int i = start; i < end; i++) {
        	GenotypingProject proj = listProj.get(i);
            String projId = Integer.toString(proj.getId());
            List<VariantSetMetadata> metadata = getMetadataList(module, projId);
            if (proj.getDescription() != null)
            {
            	VariantSetMetadata vsmd = new VariantSetMetadata();
            	vsmd.setKey("description");
            	vsmd.setValue(proj.getDescription());
            	metadata.add(vsmd);
            }
            VariantSet variantSet = VariantSet.newBuilder()
                    .setId(createId(module, projId))
                    .setReferenceSetId(module)
                    .setDatasetId(module)
                    .setName(listProj.get(i).getName())
                    .setMetadata(metadata) // get the metadata from vcf header
                    .build();
            listVariantSet.add(variantSet);
        }
        response = SearchVariantSetsResponse.newBuilder()
                .setVariantSets(listVariantSet)
                .setNextPageToken(nextPageToken)
                .build();

		return response;
    }

    @Override
    public GigwaSearchVariantsResponse searchVariants(SearchVariantsRequest svr) throws AvroRemoteException, GAException {

        GigwaSearchVariantsResponse response = null;
        // get extra info 
        GigwaSearchVariantsRequest gsvr = (GigwaSearchVariantsRequest) svr;
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(svr.getVariantSetId(), 2);
        if (info == null) {
            // wrong number of param or wrong module name 
        } else {        
            boolean doCount = false;
            boolean doSearch = false;
            boolean doBrowse = false;
            boolean getGT = gsvr.isGetGT();
            // always do count because we need it for pagination mechanism. 
            // As count are stored by query hash, it should not be a problem 
            switch (gsvr.getSearchMode()) {
                case 0:
                    doCount = true;
                    doSearch = false;
                    doBrowse = false;
                    break;
                case 1:
                    doCount = false;
                    doSearch = true;
                    doBrowse = true;
                    break;
                case 2:
                    doCount = false;
                    doSearch = false;
                    doBrowse = true;
                    break;
                case 3:
                    doCount = true;
                    doSearch = true;
                    doBrowse = true;
                    break;
            }
            String module = info[0];
            int projId = Integer.parseInt(info[1]);

            Long count = null;
            long globalCount;

            DBCursor cursor = null;
            try
            {
                if (doSearch) {
                    // create a temp collection to store the result of the request
                	count = findVariants(gsvr);
                }
                else if (doCount || doBrowse) {
                    count = countVariants(gsvr, doBrowse);
                }
                
                if (count > 0 && doBrowse)
                    cursor = getVariantCursor(gsvr);
            }
            catch (Exception ex)
            {
                LOG.error("Error searching variants", ex);
                throw new GAException(ex);
            }
            globalCount = count == null ? 0 : count;
            // get the cursor containing variant from previously created temp collection 
            // return null if cursor is empty 
            if (cursor != null && cursor.hasNext()) {
                // we need to get callSet name and position in the callSet list to get corresponding genotype
                // if we don't want to retrieve genotype, just send an empty individuals list? 
				Collection<GenotypingSample> samples;
				if (getGT) {
					samples = MgdbDao.getSamplesForProject(module, projId, gsvr.getCallSetIds().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(GigwaGa4ghServiceImpl.ID_SEPARATOR))).collect(Collectors.toList()));
				} else {
					samples = new ArrayList<>();
				}

				List<Variant> listVar = getVariantListFromDBCursor(module, Integer.parseInt(info[1]), cursor, samples, null);
				String nextPageToken = null;

                // if there is still more result after PageSize iterations
                if (globalCount > gsvr.getPageSize() * (Integer.parseInt(gsvr.getPageToken()) + 1)) {
                    nextPageToken = Integer.toString(Integer.parseInt(gsvr.getPageToken()) + 1);
                }
                response = new GigwaSearchVariantsResponse();
                response.setNextPageToken(nextPageToken);
                response.setVariants(listVar);
                if (gsvr.getSearchMode() == 3) {
                    response.setCount(count);
                }
                cursor.close();
            } else {
                response = new GigwaSearchVariantsResponse();
                response.setNextPageToken(null);
                response.setVariants(new ArrayList<>());
                response.setCount(count);
            }
        }
        return response;
    }

    @Override
    public SearchReferencesResponse searchReferences(SearchReferencesRequest srr) throws AvroRemoteException, GAException {

        SearchReferencesResponse response = null;

        // get information from id
        String[] info = GigwaSearchVariantsRequest.getInfoFromId(srr.getReferenceSetId(), 1);
        if (info == null) {

            // wrong number of param or wrong module name 
        } else {
            String module = info[0];
            GigwaSearchReferencesRequest gsr = (GigwaSearchReferencesRequest) srr;
            int projId = -1; // default : return all sequences of a module 
            String[] array = gsr.getVariantSetId().split(GigwaGa4ghServiceImpl.ID_SEPARATOR);
            if (array.length > 1) {
                projId = Integer.parseInt(gsr.getVariantSetId().split(GigwaGa4ghServiceImpl.ID_SEPARATOR)[1]);
            }

            MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
            int start;
            int end;
            int pageSize;
            int pageToken = 0;
            String nextPageToken;

            List<Reference> listReference = new ArrayList<>();
            List<String> liste = new ArrayList<>();
            Map<String, Integer> mapSeq = new TreeMap<>(new AlphaNumericComparator());

            // allow search on checksum 
            // but here only on one checksum v0.6.1 ? 
            if (gsr.getMd5checksum() != null) {

                BasicDBObject query = new BasicDBObject();
                query.put(Sequence.FIELDNAME_CHECKSUM, gsr.getMd5checksum());
                DBObject seq = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).findOne(query);
                // listSequence.add((String) seq.get("_id"));
            } else {
                Query q = new Query();
                q.fields().include(GenotypingProject.FIELDNAME_SEQUENCES);
                if (projId != -1)
                	q.addCriteria(Criteria.where("_id").is(projId));
                List<GenotypingProject> listProj = mongoTemplate.find(q, GenotypingProject.class);
                for (int i = 0; i < listProj.size(); i++) {
                    for (String seq : listProj.get(i).getSequences()) {
                        mapSeq.put(seq, i + 1);
                    }
                }
            }
            int size = mapSeq.size();
            // if page size is not specified, return all results
            if (gsr.getPageSize() != null) {
                pageSize = gsr.getPageSize();
            } else {
                pageSize = size;
            }
            if (gsr.getPageToken() != null) {
                pageToken = Integer.parseInt(gsr.getPageToken());
            }
            start = pageSize * pageToken;
            if (size - start <= pageSize) {
                end = size;
                nextPageToken = null;
            } else {
            	end = pageSize * (pageToken + 1);
                nextPageToken = Integer.toString(pageToken + 1);
            }
            ArrayList<DBObject> pipeline = new ArrayList<>();
            pipeline.add(new BasicDBObject("$match", new BasicDBObject("_id", new BasicDBObject("$in", mapSeq.keySet()))));
            Cursor sqCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).aggregate(pipeline, AggregationOptions.builder().allowDiskUse(isAggregationAllowedToUseDisk()).build());

            DBObject sequence;
            Iterator<String> iteratorName = mapSeq.keySet().iterator();
            Iterator<Integer> iteratorId = mapSeq.values().iterator();

            // create and add the corresponding Reference for each sequence 
            for (int i = start; i < end; i++) {
                long length = 0L;
                String checksum = Helper.convertToMD5("");
                String name = iteratorName.next();
                String projectId = Integer.toString(iteratorId.next());

                if (sqCursor.hasNext()) {
                    sequence = sqCursor.next();
                    length = (long) sequence.get(Sequence.FIELDNAME_LENGTH);
                    checksum = (String) sequence.get(Sequence.FIELDNAME_CHECKSUM);
                }

                String id = createId(module, projectId, name);

                // Checksum : MD5 of the upper-case sequence excluding all whitespace characters 
                // length == 0 since we don't have this information in VCF files
                Reference reference = Reference.newBuilder().setId(id)
                        .setMd5checksum(checksum)
                        .setName(name)
                        .setLength(length)
                        .setSourceAccessions(liste)
                        .build();

                listReference.add(reference);
            }

            response = SearchReferencesResponse.newBuilder()
                    .setReferences(listReference)
                    .setNextPageToken(nextPageToken)
                    .build();
        }
        return response;
    }

    /**
     * return the ID of an ontology term
     *
     * @param name
     * @return
     */
    public String getOntologyId(String name) {

        return MongoTemplateManager.getOntologyMap().get(name);
    }

    /**
     * Get annotations for a specific variant Waiting for ga4gh schema work only
     * for VCF 4.2 version using ANN notation
     *
     * @param id variant ID
     * @return snpEff annotation for this variant
     */
    public VariantAnnotation getVariantAnnotation(String id) {
        VariantAnnotation.Builder variantAnnotationBuilder = VariantAnnotation.newBuilder()
            .setVariantId(id)
            .setId(id)
            .setVariantAnnotationSetId(id.substring(0, id.lastIndexOf(ID_SEPARATOR))); // which variant annotation set?

        // get information from id
        String[] info = GigwaSearchVariantsRequest.getInfoFromId(id, 3);
        if (info == null) {
            // wrong number of param or wrong module name 
        } else {
            String module = info[0];
            String variantId = info[2];
            String[] headerField;
            String header;

            // parse annotation fields
            BasicDBObject queryVarAnn = new BasicDBObject();
            BasicDBObject varAnnField = new BasicDBObject();
            queryVarAnn.put("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, variantId);
            varAnnField.put(VariantData.FIELDNAME_KNOWN_ALLELE_LIST, 1);
            varAnnField.put(VariantData.SECTION_ADDITIONAL_INFO, 1);
            DBObject variantRunDataObj = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).findOne(queryVarAnn, varAnnField);
            DBObject variantAnnotationObj = variantRunDataObj != null ? (DBObject) variantRunDataObj.get(VariantRunData.SECTION_ADDITIONAL_INFO) : null;
            if (variantAnnotationObj != null)
            {
                String ann = (String) variantAnnotationObj.get(VcfImport.ANNOTATION_FIELDNAME_ANN);
                if (ann == null)
                	ann = (String) variantAnnotationObj.get(VcfImport.ANNOTATION_FIELDNAME_CSQ);
                boolean fAnnStyle = ann != null;
                if (!fAnnStyle)
                	ann = (String) variantAnnotationObj.get(VcfImport.ANNOTATION_FIELDNAME_EFF);
                Map<String, List<String>> additionalInfo = new HashMap<>();

                String[] tableTranscriptEffect = new String[0];
                if (ann != null)
                {
                    tableTranscriptEffect = ann.split(",");
                    List<TranscriptEffect> transcriptEffectList = new ArrayList<>();
	
	                // source version is stored in the ontology map 
	                String sourceVersion = getOntologyId(Constants.VERSION) == null ? "" : getOntologyId(Constants.VERSION);
	
	                BasicDBObject fieldHeader = new BasicDBObject(Constants.INFO_META_DATA + "." + (fAnnStyle ? VcfImport.ANNOTATION_FIELDNAME_ANN : VcfImport.ANNOTATION_FIELDNAME_EFF) + "." + Constants.DESCRIPTION, 1);
	                if (fAnnStyle)
	                	fieldHeader.put(Constants.INFO_META_DATA + "." + VcfImport.ANNOTATION_FIELDNAME_CSQ + "." + Constants.DESCRIPTION, 1);
	                
	                DBCollection vcfHeaderColl = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class));
	                BasicDBList vcfHeaderQueryOrList = new BasicDBList();
	                for (String key : fieldHeader.keySet())
	                	vcfHeaderQueryOrList.add(new BasicDBObject(key, new BasicDBObject("$exists", true)));
	
	                DBObject vcfHeaderEff = vcfHeaderColl.findOne(new BasicDBObject("$or", vcfHeaderQueryOrList), fieldHeader);
	
	                ArrayList<String> headerList = new ArrayList<>();
	                LinkedHashSet<String> usedHeaderSet = new LinkedHashSet<>();
	                if (!fAnnStyle)
	                	headerList.add("Consequence");	// EFF style annotations
	                DBObject annInfo = (DBObject) ((DBObject) vcfHeaderEff.get(Constants.INFO_META_DATA)).get(fAnnStyle ? VcfImport.ANNOTATION_FIELDNAME_ANN : VcfImport.ANNOTATION_FIELDNAME_EFF);
	                if (annInfo == null && fAnnStyle)
	                	annInfo = (DBObject) ((DBObject) vcfHeaderEff.get(Constants.INFO_META_DATA)).get(VcfImport.ANNOTATION_FIELDNAME_CSQ);
	                if (annInfo != null) {
	                    header = (String) annInfo.get(Constants.DESCRIPTION);
	                    if (header != null) {
	                        // consider using the headers for additional info keySet
	                    	String sBeforeFieldList = fAnnStyle ? ": " : " (";
	                        headerField = header.substring(header.indexOf(sBeforeFieldList) + sBeforeFieldList.length(), fAnnStyle ? header.length() : header.indexOf(")")).split("\\|");
	                        for (String head : headerField)
	                            headerList.add(head.replace("[", "").replace("]", "").trim());
	                    }
	                }
	
	                List<AnalysisResult> listAnalysisResults = new ArrayList<>();
	
	                for (int i=0; i<tableTranscriptEffect.length; i++) {
	                	ArrayList<String> values = new ArrayList<>();
	                	
	                    if (!fAnnStyle) {	// EFF style annotations
	                    	int parenthesisPos = tableTranscriptEffect[i].indexOf("(");
	                    	values.add(tableTranscriptEffect[i].substring(0, parenthesisPos));
	                    	tableTranscriptEffect[i] = tableTranscriptEffect[i].substring(parenthesisPos + 1).replace(")", "");
	                    }
	
	                    List<OntologyTerm> ontologyList = new ArrayList<>();
	
	                    String[] effectFields = tableTranscriptEffect[i].split("\\|", -1);
	                    for (int j=0; j<effectFields.length; j++)
	                    {
	                    	values.add(effectFields[j]);
	                    	if (effectFields[j].endsWith(")"))
	                    	{
	                    		String[] splitVal = effectFields[j].substring(0,  effectFields[j].length() - 1).split("\\(");
	                    		if (splitVal.length == 2)
		                    		try
		                    		{                    			
		    							AnalysisResult analysisResult = new AnalysisResult();
		    							analysisResult.setAnalysisId(headerList.get(j));
		    							analysisResult.setResult(splitVal[0]);
		    							analysisResult.setScore((int)(100 * Float.parseFloat(splitVal[1])));
		    							listAnalysisResults.add(analysisResult);
		                    		}
		                    		catch (NumberFormatException ignored)
		                    		{}
	                    	}
	                    }
	
	                    int impactIndex = headerList.indexOf(fAnnStyle ? "IMPACT" : "Effefct_Impact");
	                    if (impactIndex != -1)
		                {
		                    String[] impact = values.get(impactIndex).split("&");
		                    
		                    for (String anImpact : impact) {
		                        String ontologyId = getOntologyId(anImpact);
		                        if (ontologyId == null) {
		                            ontologyId = "";
		                        }
		                        OntologyTerm ontologyTerm = OntologyTerm.newBuilder()
		                                .setId(ontologyId)
		                                .setSourceName("sequence ontology")
		                                .setSourceVersion(sourceVersion)
		                                .setTerm(anImpact)
		                                .build();
		                        ontologyList.add(ontologyTerm);
		                    }
		                }
	
	                    HGVSAnnotation.Builder hgvsBuilder = HGVSAnnotation.newBuilder();
	                    AlleleLocation cDnaLocation = null;
	                    AlleleLocation cdsLocation = null;
	                    AlleleLocation proteinLocation = null;
	                    int nC = headerList.indexOf("HGVSc"), nP = headerList.indexOf("HGVSp"), nT = headerList.indexOf("Transcript");
	                    if ((nC != -1 && !values.get(nC).isEmpty()) || (nP != -1 && !values.get(nP).isEmpty()) || (nT != -1 && !values.get(nT).isEmpty()))
	                	{
	                		if (nC != -1)
	                			hgvsBuilder.setGenomic(values.get(nC));
	                		if (nT != -1)
	                			hgvsBuilder.setTranscript(values.get(nT));
	                		if (nP != -1)
	                			hgvsBuilder.setProtein(values.get(nP));
	                	}
	
	                    if (fAnnStyle)
	                    {
	                    	for (String positionHeader : Arrays.asList("cDNA_position", "CDS_position", "Protein_position"))
	                    	{
			                    int nPos = headerList.indexOf(positionHeader);
			                    if (nPos != -1)
			                    {
				                    String value = values.get(nPos);
			                    	if (!value.equals(""))
				                    {
				                    	AlleleLocation.Builder allLocBuilder = AlleleLocation.newBuilder();
			                    		String[] splitVals = value.split("/");
			                    		if (splitVals.length == 1 && value.contains("-"))
			                    			splitVals = value.split("-");	// sometimes used as separator
			                    		try
			                    		{
			                    			allLocBuilder.setStart(Integer.parseInt(splitVals[0]));
			                    		}
			                    		catch (NumberFormatException ignored)
			                    		{}
		
				                    	if (allLocBuilder.getStart() == 0)
				                    		continue;
				                    	
				                    	boolean fWorkingOnProtein = "Protein_position".equals(positionHeader);
		
				                    	String sRefAllele = ((List<String>) variantRunDataObj.get(VariantData.FIELDNAME_KNOWN_ALLELE_LIST)).get(0);
				                    	if (!fWorkingOnProtein)
			                    			allLocBuilder.setEnd(allLocBuilder.getStart() + sRefAllele.length() - 1);
//				                    	else
				                    		/* TODO: don't know how to calculate END field for proteins */
	
				                    	if ("cDNA_position".equals(positionHeader))
				                    		cDnaLocation = allLocBuilder.build();
				                    	else if (!fWorkingOnProtein)
				                    		cdsLocation = allLocBuilder.build();
				                    	else
				                    		proteinLocation = allLocBuilder.build();
				                    }
			                    }
	                    	}
	                    }
	
	                    TranscriptEffect transcriptEffect = TranscriptEffect.newBuilder()
	                            .setAlternateBases(values.get(0))
	                            .setId(id + ID_SEPARATOR + i)
	                            .setEffects(ontologyList)
	                            .setHgvsAnnotation(hgvsBuilder.build())
	                            .setCDNALocation(cDnaLocation)
	                            .setCDSLocation(cdsLocation)
	                            .setProteinLocation(proteinLocation)
	                            .setFeatureId(values.get(6))
	                            .setAnalysisResults(listAnalysisResults)
	                            .build();
	
	                    transcriptEffectList.add(transcriptEffect);
	                    additionalInfo.put(Constants.ANN_VALUE_LIST_PREFIX + i, values);
	                    for (int j=0; j<values.size(); j++)
	                    	if (!values.get(j).isEmpty())
	                    		usedHeaderSet.add(headerList.get(j));
	                }
	                
	                for (int i=0; i<tableTranscriptEffect.length; i++)
	                {
	                	List<String> keptValues = new ArrayList<String>(), allValues = additionalInfo.get(Constants.ANN_VALUE_LIST_PREFIX + i);
		                for (int j=0; j<allValues.size(); j++)
		                	if (usedHeaderSet.contains(headerList.get(j)))
		                		keptValues.add(allValues.get(j));
		                additionalInfo.put(Constants.ANN_VALUE_LIST_PREFIX + i, keptValues);
	                }

	                ArrayList<String> properlySortedUsedHeaderList = new ArrayList<>();
	                for (String aHeader : headerList)
	                	if (usedHeaderSet.contains(aHeader))
	                		properlySortedUsedHeaderList.add(aHeader);
	                additionalInfo.put(Constants.ANN_HEADER, new ArrayList<String>(properlySortedUsedHeaderList));
	                
	                variantAnnotationBuilder.setTranscriptEffects(transcriptEffectList);
	            }

                TreeMap<String, String> metadata = new TreeMap<>();
                for (String key : variantAnnotationObj.keySet())
                	// do not store EFF_ge / EFF_nm / EFF / ANN / CSW
                    if (!key.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE) && !key.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME) && !key.equals(VcfImport.ANNOTATION_FIELDNAME_ANN) && !key.equals(VcfImport.ANNOTATION_FIELDNAME_CSQ) && !key.equals(VcfImport.ANNOTATION_FIELDNAME_EFF) && !key.equals(""))
                        metadata.put(key, variantAnnotationObj.get(key).toString());
                additionalInfo.put(Constants.METADATA_HEADER, new ArrayList<String>(metadata.keySet()));
                additionalInfo.put(Constants.METADATA_VALUE_LIST, new ArrayList<String>(metadata.values()));
                variantAnnotationBuilder.setInfo(additionalInfo);
            }
        }
        return variantAnnotationBuilder.build();
    }

//    @Override
//    public String getReferenceBases(String seqName, int start, int end, String module) throws ObjectNotFoundException {
//
//        String sequenceBase = "";
//        int spaceCode = (int) '\n';
//        int chevCode = (int) '>';
//
//        if (end > start) {
//            Sequence seq = MongoTemplateManager.get(module).findById(seqName, Sequence.class);
//            BufferedReader br = null;
//            StringBuilder builder;
//            if ((seq != null)) {
//                try {
//                    br = new BufferedReader(new InputStreamReader(new FileInputStream(seq.getFilePath())));
//                    builder = new StringBuilder();
//                    String line;
//                    int c;
//                    int pos = 0;
//                    // skip line to reach sequence 
//                    while ((line = br.readLine()) != null && !line.startsWith(">" + seqName)) {
//
//                    }
//                    // skip char to reach start pos
//                    while ((c = br.read()) != -1 && pos < start) {
//                        pos++;
//                    }
//                    builder.append((char) c);
//
//                    while ((c = br.read()) != -1 && pos < end && c != chevCode) {
//                        if (c != spaceCode) {
//                            builder.append((char) c);
//                            pos++;
//                        }
//                    }
//                    sequenceBase = builder.toString();
//
//                } catch (IOException ex) {
//                    LOG.warn("could not open file : " + ex);
//                } finally {
//                    try {
//                        if (br != null) {
//                            br.close();
//                        }
//                    } catch (IOException ex) {
//                        LOG.warn("could not close writer : " + ex);
//                    }
//                }
//            }
//            else
//                throw new ObjectNotFoundException("No fasta for sequence " + seqName);
//        }
//        return sequenceBase;
//    }

    @Override
    public Map<String, String> getAnnotationHeaders(String module, int projId) {

        Map<String, String> annHeaders = new HashMap<>();
        BasicDBObject queryVarAnn = new BasicDBObject();
        BasicDBObject varAnnField = new BasicDBObject();
        queryVarAnn.put("_id." + DBVCFHeader.VcfHeaderId.FIELDNAME_PROJECT, projId);
        varAnnField.put(Constants.INFO_META_DATA, 1);
        varAnnField.put(Constants.INFO_FORMAT_META_DATA, 1);
        DBObject result = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).findOne(queryVarAnn, varAnnField);
        if (result != null) {
            DBObject metaDataHeader = (DBObject) result.get(Constants.INFO_META_DATA);
            for (String key : metaDataHeader.keySet()) {
                annHeaders.put(key, (String) ((DBObject) metaDataHeader.get(key)).get(Constants.DESCRIPTION));
            }
            DBObject formatHeader = (DBObject) result.get(Constants.INFO_FORMAT_META_DATA);
            for (String key : formatHeader.keySet()) {
                annHeaders.put(key, (String) ((DBObject) formatHeader.get(key)).get(Constants.DESCRIPTION));
            }
        }
        return annHeaders;
    }

    @Override
    public TreeMap<String, HashMap<String, String>> getExportFormat(String module, int projId) {
        TreeMap<String, HashMap<String, String>> exportFormats = new TreeMap<>();
        try {
            for (IExportHandler exportHandler : AbstractIndividualOrientedExportHandler.getIndividualOrientedExportHandlers().values()) {
                HashMap<String, String> info = new HashMap<>();
                info.put("desc", exportHandler.getExportFormatDescription());
                info.put("dataFileExtentions", StringUtils.join(exportHandler.getExportDataFileExtensions(), ";"));
                info.put("supportedVariantTypes", StringUtils.join(exportHandler.getSupportedVariantTypes(), ";"));
                exportFormats.put(exportHandler.getExportFormatName(), info);
            }
            for (IExportHandler exportHandler : AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().values()) {
                HashMap<String, String> info = new HashMap<>();
                info.put("desc", exportHandler.getExportFormatDescription());
                info.put("dataFileExtentions", StringUtils.join(exportHandler.getExportDataFileExtensions(), ";"));
                info.put("supportedVariantTypes", StringUtils.join(exportHandler.getSupportedVariantTypes(), ";"));
                exportFormats.put(exportHandler.getExportFormatName(), info);
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException | SecurityException ex) {
            LOG.debug("error", ex);
        }
        return exportFormats;
    }
}
