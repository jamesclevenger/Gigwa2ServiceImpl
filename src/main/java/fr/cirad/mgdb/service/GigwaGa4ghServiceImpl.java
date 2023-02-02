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

import static java.lang.Boolean.parseBoolean;

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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.HttpSession;

import org.apache.avro.AvroRemoteException;
import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;
import org.bson.Document;
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
import org.springframework.beans.factory.annotation.Value;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;
import org.springframework.security.core.Authentication;
import org.springframework.stereotype.Component;

import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.DBObject;
import com.mongodb.MongoCommandException;
import com.mongodb.client.FindIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;
import com.mongodb.client.model.Aggregates;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.exporting.individualoriented.AbstractIndividualOrientedExportHandler;
import fr.cirad.mgdb.exporting.markeroriented.AbstractMarkerOrientedExportHandler;
import fr.cirad.mgdb.importing.SequenceImport;
import fr.cirad.mgdb.importing.VcfImport;
import fr.cirad.mgdb.model.mongo.maintypes.CachedCount;
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
import fr.cirad.model.GigwaSearchCallSetsRequest;
import fr.cirad.model.GigwaSearchReferencesRequest;
import fr.cirad.model.GigwaSearchVariantsExportRequest;
import fr.cirad.model.GigwaSearchVariantsRequest;
import fr.cirad.model.GigwaSearchVariantsResponse;
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
import java.util.logging.Level;
import javax.ejb.ObjectNotFoundException;

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
public class GigwaGa4ghServiceImpl implements IGigwaService, VariantMethods, ReferenceMethods {

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

    /**
     * number of simultaneous query threads
     */
    static final int INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS = 10;
    static final private int MINIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS = 5;
    static final private int MAXIMUM_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS = 50;

    static final protected HashMap<String, String> annotationField = new HashMap<>();

    private Boolean fAllowDiskUse = null;

    @Autowired AbstractTokenManager tokenManager;

    @Autowired private AppConfig appConfig;

    private HashSet<String> hostsNotSupportingMergeOperator = new HashSet<>();

    @Autowired private MgdbDao mgdbDao;
    
    @Value("${DOCUMENT_DB_COMPAT_MODE}")
    private String documentDbCompatMode;

    public static final Integer QUERY_IDS_CHUNK_SIZE = 100000;

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
        List<String> variantTypesArray = new ArrayList<>(getProjectVariantTypes(sModule, projId));
        Collections.sort(variantTypesArray, new AlphaNumericComparator<String>());
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
        return mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(GenotypingProject.class)).distinct(GenotypingProject.FIELDNAME_ALLELE_COUNTS, projId == null ? null : new BasicDBObject("_id", projId), Integer.class).into(new ArrayList<>());
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
        List<String> indArray = null;
        try {
            indArray = new ArrayList(MgdbDao.getProjectIndividuals(sModule, project));
        } catch (ObjectNotFoundException ex) {
            java.util.logging.Logger.getLogger(GigwaGa4ghServiceImpl.class.getName()).log(Level.SEVERE, null, ex);
        }
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
    public Collection<BasicDBList> buildVariantDataQuery(GigwaSearchVariantsRequest gsvr, List<String> externallySelectedSeqs, boolean fForBrowsing) {
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
        List<String> selectedVariantIds = gsvr.getSelectedVariantIds().length() == 0 ? null : Arrays.asList(gsvr.getSelectedVariantIds().split(";"));

        Collection<BasicDBList> queries = new ArrayList<>();

        if (selectedVariantIds == null || fForBrowsing) {
            BasicDBList variantFeatureFilterList = new BasicDBList();
            /* Step to match selected variant types */
            if (selectedVariantTypes != null && selectedVariantTypes.size() > 0) {
                BasicDBList orList1 = new BasicDBList();
                BasicDBObject orSelectedVariantTypesList = new BasicDBObject();
                for (String aSelectedVariantTypes : selectedVariantTypes) {
                    BasicDBObject orClause1 = new BasicDBObject(VariantData.FIELDNAME_TYPE, aSelectedVariantTypes);
                    orList1.add(orClause1);
                    orSelectedVariantTypesList.put("$or", orList1);
                }
                variantFeatureFilterList.add(orSelectedVariantTypesList);
            }

            /* Step to match variants position range for visualization (differs from the 2 next, which is for defining the subset of data Gigwa is currently working with: it they are contradictory it still makes sense and means user is trying to view variants outside the range selected in Gigwa) */
            if (GigwaDensityRequest.class.isAssignableFrom(gsvr.getClass())) {
                variantFeatureFilterList.add(new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, ((GigwaDensityRequest) gsvr).getDisplayedSequence()));

                BasicDBObject posCrit = new BasicDBObject();
                Long min = ((GigwaDensityRequest) gsvr).getDisplayedRangeMin(), max = ((GigwaDensityRequest) gsvr).getDisplayedRangeMax();
                if (min != null && min != -1)
                    posCrit.put("$gte", min);
                if (max != null && max != -1)
                    posCrit.put("$lte", max);
                if (!posCrit.isEmpty())
                variantFeatureFilterList.add(new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, posCrit));
            }

            /* Step to match selected chromosomes */
            if (selectedSequences != null && selectedSequences.size() > 0 && selectedSequences.size() != getProjectSequences(sModule, projId).size()) {
                variantFeatureFilterList.add(new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, new BasicDBObject("$in", selectedSequences)));
            }

            /* Step to match variants that have a position included in the specified range */
            if ((gsvr.getStart() != null && gsvr.getStart() != -1) || (gsvr.getEnd() != null && gsvr.getEnd() != -1)) {
                BasicDBObject posCrit = new BasicDBObject();
                if (gsvr.getStart() != null && gsvr.getStart() != -1)
                    posCrit.put("$or", Arrays.asList(
                    	new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, new BasicDBObject("$gte", gsvr.getStart())), 
                		new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_END_SITE, new BasicDBObject("$gte", gsvr.getStart()))
                	));
                if (gsvr.getEnd() != null && gsvr.getEnd() != -1)
                    posCrit.put(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, new BasicDBObject("$lte", gsvr.getEnd()));
                variantFeatureFilterList.add(posCrit);
            }

            /* Step to match selected number of alleles */
            if (alleleCountList != null) {
                BasicDBList orList3 = new BasicDBList();
                BasicDBObject orSelectedNumberOfAllelesList = new BasicDBObject();
                for (String aSelectedNumberOfAlleles : alleleCountList) {
                    int alleleNumber = Integer.parseInt(aSelectedNumberOfAlleles);
                    orList3.add(new BasicDBObject(VariantData.FIELDNAME_KNOWN_ALLELES, new BasicDBObject("$size", alleleNumber)));
                    orSelectedNumberOfAllelesList.put("$or", orList3);
                }
                variantFeatureFilterList.add(orSelectedNumberOfAllelesList);
            }

            if (!variantFeatureFilterList.isEmpty())
                queries.add(variantFeatureFilterList);
        }
        else {    // filtering on variant IDs: we might need to split the query in order to avoid reaching a 16Mb document size
            int step = selectedVariantIds.size() / QUERY_IDS_CHUNK_SIZE;
            int r = selectedVariantIds.size() % QUERY_IDS_CHUNK_SIZE;
            if (r != 0)
                step++;

            for (int i = 0; i<step; i++) {
                List<String> subList = selectedVariantIds.subList(i*QUERY_IDS_CHUNK_SIZE, Math.min((i+1)*QUERY_IDS_CHUNK_SIZE, selectedVariantIds.size()));
                BasicDBList variantFeatureFilterList = new BasicDBList();
                variantFeatureFilterList.add(new BasicDBObject("_id", new BasicDBObject("$in", subList)));
                queries.add(variantFeatureFilterList);
            }
        }
        return queries;
    }

    private String isSearchedDatasetReasonablySized(GigwaSearchVariantsRequest gsvr) throws Exception
    {
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        int nProjectIndCount = MgdbDao.getProjectIndividuals(sModule, projId).size();
        int nGroup1IndCount = gsvr.getCallSetIds() != null && gsvr.getCallSetIds().size() != 0 ? gsvr.getCallSetIds().size() : nProjectIndCount;
        int nGroup2IndCount = gsvr.getCallSetIds2() != null && gsvr.getCallSetIds2().size() != 0 ? gsvr.getCallSetIds().size() : nProjectIndCount;

        List<Integer> groupsForWhichToFilterOnGenotypingData = GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingData(gsvr, false);
        int nIndCount = (groupsForWhichToFilterOnGenotypingData.contains(0) ? nGroup1IndCount : 0)  + (groupsForWhichToFilterOnGenotypingData.contains(1) ? nGroup2IndCount : 0);
        if (nIndCount == 0)
            return null;    // no genotyping data filtering involved

        int nMaxBillionGenotypesInvolved = 1;    // default
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
            return null;    // we can't expect user to select less than a single sequence

        int nAvgVariantsPerSeq = (int) (mongoTemplate.count(new Query(), VariantData.class) / Math.max(1, proj.getSequences().size()));
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

        boolean fGotTokenManager = tokenManager != null;    // if null, we are probably being invoked via unit-test
        String token = !fGotTokenManager ? Helper.convertToMD5(String.valueOf(System.currentTimeMillis())) /* create a mock up token */ : tokenManager.readToken(gsvr.getRequest());

        ProgressIndicator progress = ProgressIndicator.get(token);    // it may already exist (if we're being called by findVariants for example)
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

        String queryKey = getQueryKey(gsvr);
        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Long count = CachedCount.getCachedCount(mongoTemplate, queryKey, null);
        LOG.debug((count == null ? "new" : "existing") + " queryKey hash: " + queryKey);
        if (count == null)
        {
            long before = System.currentTimeMillis();
            progress.addStep("Counting matching variants");

            List<String> alleleCountList = gsvr.getAlleleCount().length() == 0 ? null : Arrays.asList(gsvr.getAlleleCount().split(";"));

            GenotypingProject genotypingProject = mongoTemplate.findById(projId, GenotypingProject.class);
            if (genotypingProject.getAlleleCounts().size() != 1 || genotypingProject.getAlleleCounts().iterator().next() != 2) {    // Project does not only have bi-allelic data: make sure we can apply MAF filter on selection
                boolean fExactlyOneNumberOfAllelesSelected = alleleCountList != null && alleleCountList.size() == 1;
                boolean fBiAllelicSelected = fExactlyOneNumberOfAllelesSelected && "2".equals(alleleCountList.get(0));
                boolean fMafRequested = (gsvr.getMaxMaf() != null && gsvr.getMaxMaf() < 50) || (gsvr.getMinMaf() != null && gsvr.getMinMaf() > 0);
                if (fMafRequested && !fBiAllelicSelected) {
                    progress.setError("MAF is only supported on biallelic data!");
                    return 0l;
                }
            }

            MongoCollection<Document> varColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
            List<Integer> filteredGroups = GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsvr, false);
            Collection<BasicDBList> variantQueryDBListColl = buildVariantDataQuery(gsvr, !fGotTokenManager ? null : getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), sModule), false);

            if (variantQueryDBListColl.isEmpty()) {
                if (filteredGroups.size() == 0 && mongoTemplate.count(new Query(), GenotypingProject.class) == 1)
                    count = mongoTemplate.count(new Query(), VariantData.class);    // no filter whatsoever
            }
            else if (filteredGroups.size() == 0) {    // filtering on variant features only: we just need a count
                count = 0l;
                count = variantQueryDBListColl.parallelStream()
                        .map(req -> varColl.countDocuments(new BasicDBObject("$and", req)))
                        .reduce(count, (accumulator, _item) -> accumulator + _item);
            }

            if (count != null)
            	mongoTemplate.save(new CachedCount(queryKey, Arrays.asList(count)));
            else
            {    // filter on genotyping data
                boolean fPreFilterOnVarColl = false, fMongoOnSameServer = MongoTemplateManager.isModuleOnLocalHost(sModule);

                //in this case, there is only one variantQueryDBList (no filtering on variant ids)
                BasicDBList variantQueryDBList = !variantQueryDBListColl.isEmpty() ? variantQueryDBListColl.iterator().next() : new BasicDBList();

                if (variantQueryDBList.size() > 0)
                {
                    Number avgObjSize = (Number) mongoTemplate.getDb().runCommand(new BasicDBObject("collStats", mongoTemplate.getCollectionName(VariantRunData.class))).get("avgObjSize");
                    if (avgObjSize.doubleValue() >= 10240)
                    {    // it may be worth pre-filtering data on variant collection because filtering speed on the run collection is affected by the document size
                        long totalCount = mongoTemplate.count(new Query(), VariantData.class), preFilterCount = varColl.countDocuments(new BasicDBObject("$and", variantQueryDBList));
                        fPreFilterOnVarColl = preFilterCount <= totalCount*(fMongoOnSameServer ? .85 : .45);    // only pre-filter if less than a given portion of the total variants are to be retained
                        if (fPreFilterOnVarColl)
                                LOG.debug("Pre-filtering data on variant collection");
                    }
                }

                GenotypingDataQueryBuilder genotypingDataQueryBuilder = new GenotypingDataQueryBuilder(gsvr, variantQueryDBList, true);
                final int nChunkCount = genotypingDataQueryBuilder.getNumberOfQueries();
                final List<Integer> shuffledChunkIndexes = genotypingDataQueryBuilder.suffleChunkOrder();

                try
                {
                    if (nChunkCount > 1)
                        LOG.debug("Query split into " + nChunkCount);

                    final Long[] partialCountArray = new Long[nChunkCount];
                    final ArrayList<Thread> threadsToWaitFor = new ArrayList<>();
                    final AtomicInteger finishedThreadCount = new AtomicInteger(0);

                    int i = -1, nNConcurrentThreads = INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS;
                    while (genotypingDataQueryBuilder.hasNext()) {
                        final List<BasicDBObject> genotypingDataPipeline = genotypingDataQueryBuilder.next();

                        final int chunkIndex = shuffledChunkIndexes.get(++i);
                        boolean fMultiProjectDB = false;

                        BasicDBObject initialMatch = (BasicDBObject) genotypingDataPipeline.get(0).get("$match");
                        if (initialMatch != null && fPreFilterOnVarColl)
                        {    // initialMatchForVariantColl will be the one applied to variants collection when pre-filtering
                            BasicDBList initialMatchForVariantColl = (BasicDBList) ((BasicDBList) initialMatch.get("$and")).clone();
                            if (initialMatchForVariantColl != null)
                            {
                                List<DBObject> toAdd = new ArrayList<>(), toRemove = new ArrayList<>();
                                for (Object filter : initialMatchForVariantColl)
                                {
                                        Object variantIdFilter = ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
                                        if (variantIdFilter != null)
                                        {
                                                toAdd.add(new BasicDBObject("_id", variantIdFilter));
                                                toRemove.add((DBObject) filter);
                                        }
                                        else if (null != ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID)) {
                                                toRemove.add((DBObject) filter);    // no project info to filter on in the variants collection
                                                fMultiProjectDB = true;
                                        }
                                }
                                initialMatchForVariantColl.addAll(toAdd);
                                initialMatchForVariantColl.removeAll(toRemove);
                            }

                            if (fMongoOnSameServer)
                            {    // always worth pre-filtering
	                            MongoCursor<Document> variantCursor = varColl.find(new BasicDBObject("$and", initialMatchForVariantColl)).projection(new BasicDBObject("_id", 1)).iterator();
	                            List<Comparable> chunkPreFilteredIDs = new ArrayList<>();
	                            while (variantCursor.hasNext())
	                            	chunkPreFilteredIDs.add((Comparable) variantCursor.next().get("_id"));
	                            if (chunkPreFilteredIDs.size() == 0)
	                            {    // no variants match indexed part of the query: skip chunk
                                    partialCountArray[chunkIndex] = 0l;
                                    // do as if one more async thread was launched so we keep better track of the progress
                                    threadsToWaitFor.add(null);
                                    finishedThreadCount.incrementAndGet();
                                    continue;
	                            }
	                            else
	                            {    // DB server is the same machine as web server: $in operator will not be expensive
                                    if (!fMultiProjectDB)    // for single project dbs, $in is equivalent to original query, otherwise only a pre-filter
                                    	genotypingDataPipeline.remove(0);
                                    genotypingDataPipeline.add(0, new BasicDBObject("$match", new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", chunkPreFilteredIDs))));
	                            }
                            }
                            else
                            {    // only try and use pre-filtering to avoid executing genotyping data queries on irrelevant chunks
                                if (varColl.countDocuments(new BasicDBObject("$and", initialMatchForVariantColl)) == 0)
                                {    // no variants match indexed part of the query: skip chunk
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

                        if (progress.isAborted())
                            return 0l;

                        final ProgressIndicator finalProgress = progress;
                        Thread queryThread = new Thread() {
                            @Override
                            public void run() {
                                    try {
                                        MongoCursor<Document> it = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(genotypingDataPipeline).allowDiskUse(isAggregationAllowedToUseDisk()).iterator();
                                        partialCountArray[chunkIndex] = it.hasNext() ? ((Number) it.next().get("count")).longValue() : 0;
                                        finalProgress.setCurrentStepProgress((short) (finishedThreadCount.incrementAndGet() * 100 / nChunkCount));
                                        genotypingDataPipeline.clear();    // release memory (VERY IMPORTANT)
                                        it.close();
                                    }
                                    catch (Throwable t) {
                                        LOG.error("Error counting variants", t);
                                        finalProgress.setError(t.getMessage());
                                        return;
                                    }
                            }
                        };

                        if (chunkIndex % nNConcurrentThreads == (nNConcurrentThreads - 1)) {
                            threadsToWaitFor.add(queryThread); // only needed to have an accurate count
                            queryThread.run();    // run synchronously

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
                            queryThread.start();    // run asynchronously for better speed
                        }
                    }

                    for (Thread t : threadsToWaitFor) // wait for all threads before moving to next phase
                        if (t != null)
                                t.join();
                    progress.setCurrentStepProgress(100);

                    count = 0l;
                    for (Long partialCount : partialCountArray)
                        count += partialCount;

                    mongoTemplate.save(new CachedCount(queryKey, Arrays.asList(partialCountArray)));
                }
                catch (InterruptedException e) {
                    LOG.debug("InterruptedException", e);
                }
            }
            LOG.info("countVariants found " + count + " results in " + (System.currentTimeMillis() - before) / 1000d + "s");
        }

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
    public MongoCollection<Document> getTemporaryVariantCollection(String sModule, String processID, boolean fEmptyItBeforeHand) {
        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        MongoCollection<Document> tmpColl = mongoTemplate.getCollection(MongoTemplateManager.TEMP_COLL_PREFIX + Helper.convertToMD5(processID));
        if (fEmptyItBeforeHand) {

//            ArrayList<StackTraceElement> keptStackTraceElements = new ArrayList<>();
//            Exception e = new Exception("Check stack trace");
//            for (StackTraceElement ste : e.getStackTrace())
//                if (ste.toString().startsWith("fr.cirad."))
//                    keptStackTraceElements.add(ste);
//            e.setStackTrace(keptStackTraceElements.toArray(new StackTraceElement[keptStackTraceElements.size()]));
//            LOG.debug("Dropping " + sModule + "." + tmpColl.getName() + " from getTemporaryVariantCollection", e);

            tmpColl.drop();
            MgdbDao.ensurePositionIndexes(mongoTemplate, Arrays.asList(tmpColl));    // make sure we have indexes defined as required in v2.4
        }
        return tmpColl;
    }

    @Override
    public long findVariants(GigwaSearchVariantsRequest gsvr) throws Exception {
        String token = tokenManager.readToken(gsvr.getRequest());

        final ProgressIndicator progress = new ProgressIndicator(token, new String[0]);
        ProgressIndicator.registerProgressIndicator(progress);
        String sizeProblemMsg = gsvr.shallApplyMatrixSizeLimit() ? isSearchedDatasetReasonablySized(gsvr) : null;
        if (sizeProblemMsg != null) {
            progress.setError(sizeProblemMsg);
            return 0;
        }

        progress.addStep("Finding matching variants");

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvr.getVariantSetId(), 2);
        String sModule = info[0];
        String queryKey = getQueryKey(gsvr);

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        String sMongoHost = MongoTemplateManager.getModuleHost(sModule);

        MongoCollection<Document> cachedCountCollection = mongoTemplate.getCollection(mongoTemplate.getCollectionName(CachedCount.class));
        //cachedCountCollection.drop();
        MongoCursor<Document> countCursor = cachedCountCollection.find(new BasicDBObject("_id", queryKey)).iterator();

        final Object[] partialCountArray = !countCursor.hasNext() ? null : ((List<Object>) countCursor.next().get(MgdbDao.FIELD_NAME_CACHED_COUNT_VALUE)).toArray();
        final LinkedHashMap<Integer, Long> partialCountMap = new LinkedHashMap<>(); // progress display will be more accurate if we skip empty chunks
        long nTotalCount = 0;
        if (partialCountArray != null) {
            for (int i=0; i<partialCountArray.length; i++) {
                long n = (long) partialCountArray[i];
                if (n != 0) {
                    partialCountMap.put(i, n);
                    nTotalCount += n;
                }
            }
            if (nTotalCount == 0) {
                progress.markAsComplete();
                return 0;
            }
        }
        LOG.debug((partialCountArray == null ? "new" : "existing") + " queryKey hash: " + queryKey);

        List<Integer> filteredGroups = GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsvr, false);
        boolean fNeedToCreateTempColl = !filteredGroups.isEmpty() || (gsvr.getSelectedVariantIds() != null && !gsvr.getSelectedVariantIds().isEmpty());        

        long before = System.currentTimeMillis();

        if (fNeedToCreateTempColl) {
            final MongoCollection<Document> tmpVarColl = getTemporaryVariantCollection(sModule, progress.getProcessId(), true);
            Collection<BasicDBList> variantQueryDBListColl = buildVariantDataQuery(gsvr, getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), sModule), false);

	        if (filteredGroups.isEmpty()) {	// filtering by variant IDs
	            MongoCollection<Document> varColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
	            final AtomicInteger nProcessedChunkCount = new AtomicInteger(0);
	
	            variantQueryDBListColl.parallelStream().forEach(req -> {
	                varColl.aggregate(Arrays.asList(
	                    new BasicDBObject("$match", new BasicDBObject("$and", req)),
	                    new BasicDBObject("$merge", new BasicDBObject("into", tmpVarColl.getNamespace().getCollectionName()).append("whenMatched", "fail"))
	                )).toCollection();
	
	                progress.setCurrentStepProgress(nProcessedChunkCount.incrementAndGet() * 100 / variantQueryDBListColl.size());
	            });
	        }
	        else {   // filter on genotyping data
	            final ArrayList<Thread> threadsToWaitFor = new ArrayList<>();
	            final AtomicInteger finishedThreadCount = new AtomicInteger(0);
	
	            //in this case, there is only one variantQueryDBList (no filtering on variant ids)
	            BasicDBList variantQueryDBList = !variantQueryDBListColl.isEmpty() ? variantQueryDBListColl.iterator().next() : new BasicDBList();
	
	            final GenotypingDataQueryBuilder genotypingDataQueryBuilder = new GenotypingDataQueryBuilder(gsvr, variantQueryDBList, false);
	
	            try {
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
	
	                MongoCollection<Document> varColl = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class));
	
	                boolean fPreFilterOnVarColl = false, fMongoOnSameServer = MongoTemplateManager.isModuleOnLocalHost(sModule);
	                if (variantQueryDBList.size() > 0) {
	                        Number avgObjSize = (Number) mongoTemplate.getDb().runCommand(new BasicDBObject("collStats", mongoTemplate.getCollectionName(VariantRunData.class))).get("avgObjSize");
	                    if (avgObjSize.doubleValue() >= 10240) {   // it may be worth pre-filtering data on variant collection because filtering speed on the run collection is affected by the document size
	                        long totalCount = mongoTemplate.count(new Query(), VariantData.class), preFilterCount = varColl.countDocuments(new BasicDBObject("$and", variantQueryDBList));
	                        fPreFilterOnVarColl = preFilterCount <= totalCount*(fMongoOnSameServer ? .85 : .45);    // only pre-filter if less than a given portion of the total variants are to be retained
	                        if (fPreFilterOnVarColl)
	                            LOG.debug("Pre-filtering data on variant collection");
	                    }
	                }
	
	                if (nChunkCount > 1)
	                    LOG.debug("Query split into " + nChunkCount);
	
	                int i = -1, nNConcurrentThreads = INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS;
	                final MongoCollection<Document> vrdColl = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class));
	                final HashMap<Integer, BasicDBList> rangesToCount = new HashMap<>();
	                while (genotypingDataQueryBuilder.hasNext()) {
	                    List<BasicDBObject> genotypingDataPipeline = genotypingDataQueryBuilder.next();
	                    if (progress.isAborted() || progress.getError() != null)
	                        return 0;
	
	                    final int chunkIndex = shuffledChunkIndexes.get(++i);
	                    if (partialCountMap.size() > 0 && !partialCountMap.containsKey(chunkIndex))
	                        continue;
	
	                    boolean fMultiProjectDB = false;
	
	                    BasicDBObject initialMatch = (BasicDBObject) genotypingDataPipeline.get(0).get("$match");
	                    BasicDBList initialMatchForVariantColl = (BasicDBList) ((BasicDBList) initialMatch.get("$and")).clone();
	                    rangesToCount.put(chunkIndex, initialMatchForVariantColl);
	                    List<DBObject> toAdd = new ArrayList<>(), toRemove = new ArrayList<>();
	                    for (Object filter : initialMatchForVariantColl)
	                    {
	                        Object variantIdFilter = ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID);
	                        if (variantIdFilter != null)
	                        {
	                            toAdd.add(new BasicDBObject("_id", variantIdFilter));
	                            toRemove.add((DBObject) filter);
	                        }
	                        else if (null != ((DBObject) filter).get("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID)) {
	                            toRemove.add((DBObject) filter);    // no project info to filter on in the variants collection
	                            fMultiProjectDB = true;
	                        }
	                    }
	                    initialMatchForVariantColl.addAll(toAdd);
	                    initialMatchForVariantColl.removeAll(toRemove);
	
	                    if (fPreFilterOnVarColl)
	                    {
	                        if (fMongoOnSameServer)
	                        {   // always worth pre-filtering
	                            MongoCursor<Document> variantCursor = varColl.find(new BasicDBObject("$and", initialMatchForVariantColl)).projection(new BasicDBObject("_id", 1)).iterator();
	                            List<Comparable> chunkPreFilteredIDs = new ArrayList<>();
	                            while (variantCursor.hasNext())
	                                chunkPreFilteredIDs.add((Comparable) variantCursor.next().get("_id"));
	                            if (chunkPreFilteredIDs.size() == 0)
	                            {   // no variants match indexed part of the query: skip chunk
	                                if (partialCountArrayToFill != null)
	                                    partialCountArrayToFill[chunkIndex] = 0l;
	                                // do as if one more async thread was launched so we keep better track of the progress
	                                threadsToWaitFor.add(null);
	                                finishedThreadCount.incrementAndGet();
	                                continue;
	                            }
	                            else
	                            {   // DB server is the same machine as web server: $in operator will not be expensive
	                                if (!fMultiProjectDB)   // for single project dbs, $in is equivalent to original query, otherwise only a pre-filter
	                                    genotypingDataPipeline.remove(0);
	                                genotypingDataPipeline.add(0, new BasicDBObject("$match", new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", chunkPreFilteredIDs))));
	                            }
	                        }
	                        else
	                        {   // only try and use pre-filtering to avoid executing genotyping data queries on irrelevant chunks
	                            if (varColl.countDocuments(new BasicDBObject("$and", initialMatchForVariantColl)) == 0)
	                            {   // no variants match indexed part of the query: skip chunk
	                                if (partialCountArrayToFill != null)
	                                    partialCountArrayToFill[chunkIndex] = 0l;
	                                // do as if one more async thread was launched so we keep better track of the progress
	                                threadsToWaitFor.add(null);
	                                finishedThreadCount.incrementAndGet();
	                                continue;
	                            }
	                        }
	                    }
	
	                    if (partialCountArray != null)
	                        genotypingDataPipeline.add(new BasicDBObject("$limit", partialCountArray[chunkIndex]));
	                    genotypingDataPipeline.add(new BasicDBObject("$project", new BasicDBObject(VariantData.FIELDNAME_KNOWN_ALLELES, 1).append(VariantData.FIELDNAME_REFERENCE_POSITION, 1).append(VariantData.FIELDNAME_TYPE, 1)));
	
	                        Thread queryThread = new Thread() {
	                            @Override
	                            public void run() {
	                                    boolean fMergeFailedOnThisChunk = false;
	                                    if (!hostsNotSupportingMergeOperator.contains(sMongoHost))
	                                        try {
	                                            genotypingDataPipeline.add(new BasicDBObject("$merge", new BasicDBObject("into", tmpVarColl.getNamespace().getCollectionName()).append("whenMatched", "fail" /* important (fastest option)*/)));
	//                                            System.out.println(genotypingDataPipeline);
	                                            vrdColl.aggregate(genotypingDataPipeline).allowDiskUse(isAggregationAllowedToUseDisk()).toCollection();
	                                        }
	                                        catch (Throwable t) {
	                                            if (t instanceof MongoCommandException && t.getMessage().contains("$merge")) {
	                                                    hostsNotSupportingMergeOperator.add(sMongoHost);
	                                                    fMergeFailedOnThisChunk = true;
	                                                    LOG.warn("Disabling use of $merge in creating temporary collections on host " + sMongoHost + " (operator not supported by MongoDB server version)");
	                                            }
	                                            else {
	                                                    LOG.error("Error searching variants", t);
	                                                    progress.setError(t.getMessage());
	                                                    return;
	                                            }
	                                        }
	                                if (hostsNotSupportingMergeOperator.contains(sMongoHost))
	                                    try {
	                                            if (fMergeFailedOnThisChunk)
	                                                    genotypingDataPipeline.remove(genotypingDataPipeline.size() - 1);    // remove the $merge step we added
	                                                MongoCursor<Document> genotypingDataCursor = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(genotypingDataPipeline).allowDiskUse(isAggregationAllowedToUseDisk()).iterator();
	                                                final ArrayList<Document> variantsThatPassedRunFilterForThisChunk = new ArrayList<>();
	                                                while (genotypingDataCursor.hasNext())
	                                                    variantsThatPassedRunFilterForThisChunk.add(genotypingDataCursor.next());
	
	                                                if (partialCountArrayToFill != null)
	                                                    partialCountArrayToFill[chunkIndex] = (long) variantsThatPassedRunFilterForThisChunk.size();
	                                                if (variantsThatPassedRunFilterForThisChunk.size() > 0)
	                                                    tmpVarColl.insertMany(variantsThatPassedRunFilterForThisChunk);
	
	                                                genotypingDataCursor.close();
	                                    }
	                                    catch (Exception e) {
	                                        LOG.error("Error searching variants", e);
	                                        progress.setError(e.getMessage());
	                                    }
	                                    genotypingDataPipeline.clear();    // release memory
	                                        finishedThreadCount.incrementAndGet();
	                            }
	                        };
	
	                        if (chunkIndex % nNConcurrentThreads == (nNConcurrentThreads - 1)) {
	                            threadsToWaitFor.add(queryThread); // we only need to have an accurate count
	                            queryThread.run();  // run synchronously
	
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
	                            queryThread.start();    // run asynchronously for better speed
	                        }
	                        progress.setCurrentStepProgress((short) (i * 100 / nChunkCount));
	                    }
	
	                    for (Thread t : threadsToWaitFor) // wait for all threads before moving to next phase
	                        if (t != null)
	                            t.join();
	
	                    if (progress.getError() == null) {
	                        progress.setCurrentStepProgress(100);
	
	                        if (partialCountArrayToFill != null) {    // we don't have a count cache for this query: let's create it
	                            if (!hostsNotSupportingMergeOperator.contains(sMongoHost)) {
	                            threadsToWaitFor.clear();
	                                for (Integer j : rangesToCount.keySet()) {
	                                    Thread countThread = new Thread() {
	                                        public void run() {
	                                            partialCountArrayToFill[j] = tmpVarColl.countDocuments(new BasicDBObject("$and", rangesToCount.get(j)));
	                                        }
	                                    };
	                                    threadsToWaitFor.add(countThread);
	
	                                    if (j % nNConcurrentThreads*2 == 0 || j == rangesToCount.size() - 1) {
	                                        for (Thread t : threadsToWaitFor)
	                                            t.start();
	
	                                        for (Thread t : threadsToWaitFor) // wait for all threads before moving on
	                                            t.join();
	
	                                        threadsToWaitFor.clear();
	                                    }
	                                }
	                            }
	                            mongoTemplate.save(new CachedCount(queryKey, Arrays.asList(partialCountArrayToFill)));
	                        }
	                    }
	                }
	                catch (InterruptedException e) {
	                    LOG.debug("InterruptedException : " + e);
	                    // throw e;
	                }
	        	}
        }
        if (partialCountArray == null)
            nTotalCount = countVariants(gsvr, true);
        else
            LOG.info("findVariants found " + nTotalCount + " results in " + (System.currentTimeMillis() - before) / 1000d + "s");

        if (progress.isAborted() || progress.getError() != null)
            return 0;

        progress.markAsComplete();
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
                } catch (IOException e) {
                    LOG.error("Unable to cleanup expired export files", e);
                }
            }
        }.start();

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsver.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        long before = System.currentTimeMillis();
        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        int nGroupsToFilterGenotypingDataOn = GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gsver, true).size();

        Set<String> selectedIndividualList1 = gsver.getCallSetIds().size() == 0 ? MgdbDao.getProjectIndividuals(sModule, projId) /* no selection means all selected */ : gsver.getCallSetIds().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(IGigwaService.ID_SEPARATOR))).collect(Collectors.toSet());
        Set<String> selectedIndividualList2 = gsver.getCallSetIds2().size() == 0 && nGroupsToFilterGenotypingDataOn > 1 ? MgdbDao.getProjectIndividuals(sModule, projId) /* no selection means all selected */ : gsver.getCallSetIds2().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(IGigwaService.ID_SEPARATOR))).collect(Collectors.toSet());
        Collection<String> individualsToExport = gsver.getExportedIndividuals().size() > 0 ? gsver.getExportedIndividuals() : MgdbDao.getProjectIndividuals(sModule, projId);

        long count = countVariants(gsver, true);
        MongoCollection<Document> tmpVarColl = getTemporaryVariantCollection(sModule, token, false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        Collection<BasicDBList> variantQueryDBListColl = buildVariantDataQuery(gsver, getSequenceIDsBeingFilteredOn(gsver.getRequest().getSession(), sModule), true);
        final BasicDBList variantQueryDBList = variantQueryDBListColl.size() == 1 ? variantQueryDBListColl.iterator().next() : new BasicDBList();

        if (nGroupsToFilterGenotypingDataOn > 0 && nTempVarCount == 0)
        {
            progress.setError(Constants.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return;
        }

        Authentication auth = tokenManager.getAuthenticationFromToken(tokenManager.readToken(gsver.getRequest()));
        String mainVarCollName = mongoTemplate.getCollectionName(VariantData.class), usedVarCollName = nTempVarCount == 0 ? mainVarCollName : tmpVarColl.getNamespace().getCollectionName();
        MongoCollection<Document> usedVarColl = mongoTemplate.getCollection(usedVarCollName);
        Document variantQuery = nTempVarCount == 0 && !variantQueryDBList.isEmpty() ? new Document("$and", variantQueryDBList) : new Document();
        if (gsver.shallApplyMatrixSizeLimit())
        {    // make sure the matrix is not too big
            int nMaxBillionGenotypesInvolved = 1;    // default
            try {
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

            BigInteger matrixSize = BigInteger.valueOf(usedVarColl.countDocuments(variantQuery)).multiply(BigInteger.valueOf(individualsToExport.size()));
            BigInteger maxAllowedSize = BigInteger.valueOf(1000000000).multiply(BigInteger.valueOf(nMaxBillionGenotypesInvolved));

            if (matrixSize.divide(maxAllowedSize).intValue() >= 1)
            {
                progress.setError("You may only export up to " + nMaxBillionGenotypesInvolved + " billion genotypes. The current selection contains " + BigDecimal.valueOf(matrixSize.longValue()).divide(BigDecimal.valueOf(1000000000)).setScale(2, BigDecimal.ROUND_HALF_UP) + " billion genotypes.");
                return;
            }
        }

        progress.addStep("Identifying matching variants");

        OutputStream os = null;

        try
        {
        	String username = AbstractTokenManager.getUserNameFromAuthentication(tokenManager.getAuthenticationFromToken(token));

            AbstractIndividualOrientedExportHandler individualOrientedExportHandler = AbstractIndividualOrientedExportHandler.getIndividualOrientedExportHandlers().get(gsver.getExportFormat());
            AbstractMarkerOrientedExportHandler markerOrientedExportHandler = AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().get(gsver.getExportFormat());

            String filename = sModule + "__project" + projId + "__" + new SimpleDateFormat("yyyy-MM-dd").format(new Date()) + "__" + count + "variants__" + gsver.getExportFormat() + "." + (individualOrientedExportHandler != null ? individualOrientedExportHandler : markerOrientedExportHandler).getExportArchiveExtension();

            LOG.info((gsver.isKeepExportOnServer() ? "On-server" : "Direct-download") + " export requested: " + processId);
            if (gsver.isKeepExportOnServer()) {
                String relativeOutputFolder = FRONTEND_URL + File.separator + TMP_OUTPUT_FOLDER + File.separator + username + File.separator + Helper.convertToMD5(processId) + File.separator;
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
                    Thread exportThread = new IExportHandler.SessionAttributeAwareExportThread(gsver.getRequest().getSession()) {
                        public void run() {
                            try {
                                progress.addStep("Reading and re-organizing genotypes"); // initial step will consist in organizing genotypes by individual rather than by marker
                                progress.moveToNextStep();    // done with identifying variants
                                File[] exportFiles = individualOrientedExportHandler.createExportFiles(sModule, nTempVarCount == 0 ? null : usedVarCollName, variantQuery, count, selectedIndividualList1, selectedIndividualList2, processId, gsver.getAnnotationFieldThresholds(), gsver.getAnnotationFieldThresholds2(), samplesToExport, progress);

                                for (String step : individualOrientedExportHandler.getStepList())
                                    progress.addStep(step);
                                progress.moveToNextStep();
                                individualOrientedExportHandler.exportData(finalOS, sModule, AbstractTokenManager.getUserNameFromAuthentication(auth), exportFiles, true, progress, usedVarCollName, variantQuery, count, null, gsver.getMetadataFields(), readyToExportFiles);
                                if (!progress.isAborted()) {
                                    LOG.info("exportVariants (" + gsver.getExportFormat() + ") took " + (System.currentTimeMillis() - before) / 1000d + "s to process " + count + " variants and " + individualsToExport.size() + " individuals");
                                    progress.markAsComplete();
                                }
                            }
                            catch (Exception e) {
                                LOG.error("Error exporting data", e);
                                progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
                            }
                            finally {
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
                progress.moveToNextStep();    // done with identifying variants

                String contentType = markerOrientedExportHandler.getExportContentType();
                if (contentType != null && contentType.trim().length() > 0)
                    response.setContentType(contentType);

                Thread exportThread = new IExportHandler.SessionAttributeAwareExportThread(gsver.getRequest().getSession()) {
                    public void run() {
                        try {
                            markerOrientedExportHandler.exportData(finalOS, sModule, AbstractTokenManager.getUserNameFromAuthentication(auth), selectedIndividualList1, selectedIndividualList2, progress, nTempVarCount == 0 ? null : usedVarColl.getNamespace().getCollectionName(), variantQuery, count, null, gsver.getAnnotationFieldThresholds(), gsver.getAnnotationFieldThresholds2(), samplesToExport, gsver.getMetadataFields(), readyToExportFiles);
                            if (!progress.isAborted() && progress.getError() == null) {
                                LOG.info("exportVariants (" + gsver.getExportFormat() + ") took " + (System.currentTimeMillis() - before) / 1000d + "s to process " + count + " variants and " + individualsToExport.size() + " individuals");
                                progress.markAsComplete();
                            }
                        }
                        catch (Exception e) {
                            LOG.error("Error exporting data", e);
                            progress.setError("Error exporting data: " + e.getClass().getSimpleName() + (e.getMessage() != null ? " - " + e.getMessage() : ""));
                        }
                        finally {
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
                sc.nextLine();    // skip queryKey line
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
            return;    // working around some random bug
        }
        long nowMillis = new Date().getTime();
        File filterOutputLocation = new File(request.getSession().getServletContext().getRealPath(FRONTEND_URL + File.separator + TMP_OUTPUT_FOLDER));
        if (filterOutputLocation.exists() && filterOutputLocation.isDirectory()) {
            for (File f : filterOutputLocation.listFiles()) {
                if (f.isDirectory() && nowMillis - f.lastModified() > EXPORT_EXPIRATION_DELAY_MILLIS) {
                    FileUtils.deleteDirectory(f);    // it is an expired job-output-folder
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
        String collName = getTemporaryVariantCollection(sModule, processID, false).getNamespace().getCollectionName();
        MongoTemplateManager.get(sModule).dropCollection(collName);
        LOG.debug("Dropped collection " + sModule + "." + collName);
    }

    @Override
    public Collection<String> distinctSequencesInSelection(HttpServletRequest request, String sModule, int projId, String processID) {
        String sShortProcessID = processID/*.substring(1 + processID.indexOf('|'))*/;
        MongoCollection<Document> tmpVarColl = getTemporaryVariantCollection(sModule, sShortProcessID, false);
        if (tmpVarColl.estimatedDocumentCount() == 0) {
            return listSequences(request, sModule, projId);    // working on full dataset
        }
        List<String> distinctSequences = tmpVarColl.distinct(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, String.class).into(new ArrayList<>());
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
                        + gsvr.getSelectedVariantIds() + ":"

                        + gsvr.getCallSetIds() + ":"
                        + gsvr.getAnnotationFieldThresholds() + ":"
                        + gsvr.getGtPattern() + ":"
                        + gsvr.getMostSameRatio() + ":"
                        + gsvr.getMinMissingData() + ":"
                        + gsvr.getMaxMissingData() + ":"
                        + gsvr.getMinMaf() + ":"
                        + gsvr.getMaxMaf() + ":"

                        + gsvr.getCallSetIds2() + ":"
                        + gsvr.getAnnotationFieldThresholds2() + ":"
                        + gsvr.getGtPattern2() + ":"
                        + gsvr.getMostSameRatio2() + ":"
                        + gsvr.getMinMissingData2() + ":"
                        + gsvr.getMaxMissingData2() + ":"
                        + gsvr.getMinMaf2() + ":"
                        + gsvr.getMaxMaf2() + ":"

                        + gsvr.isDiscriminate() + ":"
                        + gsvr.getVariantEffect();
//        System.out.println(queryKey);
            return Helper.convertToMD5(queryKey);
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
        Document sequence;

        // no need to filter on projId since we are searching on sequenceId
        ArrayList<BasicDBObject> aggregationParam = new ArrayList<>();
        aggregationParam.add(new BasicDBObject("$match", new BasicDBObject("_id", new BasicDBObject("$in", sequenceList))));

        MongoCursor<Document> genotypingDataCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(aggregationParam).allowDiskUse(isAggregationAllowedToUseDisk()).iterator();

        if (genotypingDataCursor != null && genotypingDataCursor.hasNext()) {

            while (genotypingDataCursor.hasNext()) {

                sequence = genotypingDataCursor.next();

                // empty the table
                info.clear();

                // no need for sequence?
                info.put("sequence", sequence.get(Sequence.FIELDNAME_SEQUENCE));

                info.put("length", sequence.get(Sequence.FIELDNAME_LENGTH));
                info.put("checksum", sequence.get(Sequence.FIELDNAME_CHECKSUM));

                listSeqInfo.put((String) sequence.get("_id"), info);

            }
        }
        return listSeqInfo;
    }

    /**
     * get a list of variant in ga4gh format from a MongoCursor<Document>
     *
     * @param module
     * @param projId
     * @param cursors
     * @param samples
     * @return List<Variant>
     * @throws AvroRemoteException
     */
    public List<Variant> getVariantListFromDBCursor(String module, int projId, MongoCursor<Document> cursor, Collection<GenotypingSample> samples) throws AvroRemoteException
    {
//        long before = System.currentTimeMillis();
        LinkedHashMap<Comparable, Variant> varMap = new LinkedHashMap<>();

        // parse the cursor to create all GAVariant
        while (cursor.hasNext()) {
            Document obj = cursor.next();
            // save the Id of each variant in the cursor
            String id = (String) obj.get("_id");
            List<String> knownAlleles = ((List<String>) obj.get(VariantData.FIELDNAME_KNOWN_ALLELES));

            Variant.Builder variantBuilder = Variant.newBuilder()
                    .setId(createId(module, projId, id.toString()))
                    .setVariantSetId(createId(module, projId));

            Document rp = ((Document) obj.get(VariantData.FIELDNAME_REFERENCE_POSITION));
            if (rp == null)
                variantBuilder.setReferenceName("").setStart(-1).setEnd(-1);
            else {
                String chr = (String) rp.get(ReferencePosition.FIELDNAME_SEQUENCE);
                    variantBuilder.setReferenceName(chr != null ? chr : "");

                    Long start = (Long) rp.get(ReferencePosition.FIELDNAME_START_SITE);
                variantBuilder.setStart(start != null ? start : -1);

                    Long end = (Long) rp.get(ReferencePosition.FIELDNAME_END_SITE);
                    if (end == null && start != null)
                        end = start;
                variantBuilder.setEnd(end != null ? end : -1);
            }
            if (knownAlleles.size() == 0)
                throw new AvroRemoteException("Variant " + id.toString() + " has no known alleles!");

            variantBuilder.setReferenceBases(knownAlleles.get(0)); // reference is the first one in VCF files
            variantBuilder.setAlternateBases(knownAlleles.subList(1, knownAlleles.size()));
            
            // add the annotation map to the variant
            Map<String, List<String>> annotations = new HashMap<>();
            List<String> infoType = new ArrayList<>();
            infoType.add((String) obj.get(VariantData.FIELDNAME_TYPE));
            annotations.put("type", infoType);
            variantBuilder.setInfo(annotations);
            
            varMap.put(id, variantBuilder.build());
        }

        // get the VariantRunData containing annotations
        ArrayList<BasicDBObject> pipeline = new ArrayList<>();
        // wanted fields
        BasicDBObject fields = new BasicDBObject();
        fields.put(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_GENE, 1);
        fields.put(VariantRunData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, 1);

        // get the genotype for wanted individuals/callSet only
        final Map<Integer, String> sampleIdToIndividualMap = new HashMap<>();
        for (GenotypingSample sample : samples){
            fields.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES + "." + sample.getId(), 1);
            sampleIdToIndividualMap.put(sample.getId(), sample.getIndividual());
        }

        BasicDBList matchAndList = new BasicDBList();
        matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, new BasicDBObject("$in", varMap.keySet())));
        matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_PROJECT_ID, projId));
        if (!samples.isEmpty())
            matchAndList.add(new BasicDBObject("_id." + VariantRunDataId.FIELDNAME_RUNNAME, new BasicDBObject("$in", samples.stream().map(sp -> sp.getRun()).distinct().collect(Collectors.toList()))));
        pipeline.add(new BasicDBObject("$match", new BasicDBObject("$and", matchAndList)));
        pipeline.add(new BasicDBObject("$project", fields));
        if (samples.isEmpty())    // if no genotypes are expected back then we assume we're building the result table (thus we need to include variant name & effect when available in one of then runs)
            pipeline.add(new BasicDBObject("$sort", new BasicDBObject(AbstractVariantData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, -1)));  // if some VariantRunData records have gene info they will appear first, which will make that info available for building the result table

        HashSet<String> variantsForWhichAnnotationWasRetrieved = new HashSet<>();

        MongoCursor<Document> genotypingDataCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).aggregate(pipeline).allowDiskUse(isAggregationAllowedToUseDisk()).iterator();
        while (genotypingDataCursor.hasNext())
        {
            Document variantObj = genotypingDataCursor.next();
            String varId = (String) Helper.readPossiblyNestedField(variantObj, "_id." + VariantRunDataId.FIELDNAME_VARIANT_ID, "; ");
            Variant var = varMap.get(varId);
            if (var == null /* should not happen! */|| variantsForWhichAnnotationWasRetrieved.contains(varId))
                continue;

            variantsForWhichAnnotationWasRetrieved.add(varId);
            TreeSet<Call> calls = new TreeSet(new AlphaNumericComparator<Call>());    // for automatic sorting

            // for each annotation field
            for (String key : variantObj.keySet()) {
                switch (key) {
                    // this goes in Call  || should not be called if sp field is not present
                    case VariantRunData.FIELDNAME_SAMPLEGENOTYPES:
                        // get genotype map
                        Map<String, Object> callMap = (Map<String, Object>) variantObj.get(key);

                        // for each individual/CallSet
                        for (Integer sampleId : sampleIdToIndividualMap.keySet()) {
                            Document callObj = (Document) callMap.get("" + sampleId);
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
                                var.getInfo().put(subKey, listGene);
                            } else if (subKey.equals(VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME)) {
                                List<String> listEffect = (List<String>) additionalInfos.get(subKey);
                                var.getInfo().put(subKey, listEffect);
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
        }

//        LOG.debug("getVariantListFromDBCursor took " + (System.currentTimeMillis() - before) / 1000f + "s for " + varMap.size() + " variants and " + samples.size() + " samples");
        return new ArrayList<Variant>(varMap.values());
    }

    /**
     * create composite ID from a list of params
     *
     * @param params
     * @return String the id
     */
    static public String createId(Comparable... params) {
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
        MongoCursor<Document> cursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).find(whereQuery).iterator();

        if (cursor != null && cursor.hasNext()) {

            while (cursor.hasNext()) {

                // get the vcf header of each run of a project
                vcfHeader = DBVCFHeader.fromDocument(cursor.next());

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
    public VariantSet getVariantSet(String id) throws AvroRemoteException
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
                    vsmd.setKey(Constants.DESCRIPTION);
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
    public Variant getVariant(String id) throws AvroRemoteException {
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
            String name = info[2];
            String run = null;
            if (info.length == 4)
                run = info[3];

            MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
            MongoCursor<Document> cursor = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantData.class)).find(new BasicDBObject("_id", name)).iterator();

            if (cursor != null && cursor.hasNext()) {
                List<Criteria> sampleQueryCriteria = new ArrayList<>();
                if (!listInd.isEmpty())
                    sampleQueryCriteria.add(Criteria.where(GenotypingSample.FIELDNAME_INDIVIDUAL).in(listInd));
                if (run != null)
                    sampleQueryCriteria.add(Criteria.where(GenotypingSample.FIELDNAME_RUN).is(run));
                Collection<GenotypingSample> samples = mongoTemplate.find(sampleQueryCriteria.isEmpty() ? new Query() : new Query(new Criteria().andOperator(sampleQueryCriteria.toArray(new Criteria[sampleQueryCriteria.size()]))), GenotypingSample.class);
                variant = getVariantListFromDBCursor(module, Integer.parseInt(info[1]), cursor, samples).get(0);
                cursor.close();
            }
        }
        return variant;
    }

    @Override
    public CallSet getCallSet(String id) throws AvroRemoteException {
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

            try {
                // check if the callSet is in the list
                if (MgdbDao.getProjectIndividuals(module, projId).contains(name))
                    callSet = CallSet.newBuilder().setId(id).setName(name).setVariantSetIds(listVariantSetId).setSampleId(null).build();
            } catch (ObjectNotFoundException ex) {
                java.util.logging.Logger.getLogger(GigwaGa4ghServiceImpl.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return callSet;
    }

    @Override
    public ReferenceSet getReferenceSet(String id) throws AvroRemoteException {
        ReferenceSet referenceSet = null;

        MongoTemplate mongoTemplate = MongoTemplateManager.get(id);
        if (mongoTemplate == null) {

        } else {
            List<String> list = new ArrayList<>();
            // get the all references of the reference Set/module
            MongoCursor<Document> genotypingDataCursor = MongoTemplateManager.get(id).getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).find().iterator();
            Document seq;

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
                    .setDescription((taxoDesc.isEmpty() ? "" : (taxoDesc + " ; ")) + mongoTemplate.getCollection(mongoTemplate.getCollectionName(GenotypingProject.class)).distinct(GenotypingProject.FIELDNAME_SEQUENCES, String.class).into(new ArrayList<>()).size() + " references ; " + mongoTemplate.count(new Query(), VariantData.class) + " markers")
                    .build();
        }
        return referenceSet;
    }

    @Override
    public Reference getReference(String id) throws AvroRemoteException {

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
    public ListReferenceBasesResponse getReferenceBases(String id, ListReferenceBasesRequest lrbr) throws AvroRemoteException
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

    @Override
    public SearchCallSetsResponse searchCallSets(SearchCallSetsRequest scsr) throws AvroRemoteException {
        // get information from id
        String[] info = GigwaSearchVariantsRequest.getInfoFromId(scsr.getVariantSetId(), 2);
        if (info == null)
            return null;

        GigwaSearchCallSetsRequest gscsr = (GigwaSearchCallSetsRequest) scsr;
        Authentication auth = tokenManager.getAuthenticationFromToken(tokenManager.readToken(gscsr.getRequest()));

        CallSet callSet;
        int start;
        int end;
        int pageSize;
        int pageToken = 0;
        String nextPageToken;

        String module = info[0];
        LinkedHashMap<String, Individual> indMap = mgdbDao.loadIndividualsWithAllMetadata(module, AbstractTokenManager.getUserNameFromAuthentication(auth), Arrays.asList(Integer.parseInt(info[1])), null);

        List<CallSet> listCallSet = new ArrayList<>();
        int size = indMap.size();
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
        List<String> indList = new ArrayList() {{ addAll(indMap.keySet()); }};
        for (int i = start; i < end; i++) {
            final Individual ind = indMap.get(indList.get(i));
            CallSet.Builder csb = CallSet.newBuilder().setId(createId(module, info[1], ind.getId())).setName(ind.getId()).setVariantSetIds(Arrays.asList(scsr.getVariantSetId())).setSampleId(createId(module, info[1], ind.getId(), ind.getId()));

            if (!ind.getAdditionalInfo().isEmpty()) {
                            Map<String, String> addInfoMap = new HashMap();
                            for (String key:ind.getAdditionalInfo().keySet()) {
                                Object value = ind.getAdditionalInfo().get(key);
                                if (value instanceof String) {
                                    int spaces = ((String) value).length() - ((String) value).replaceAll(" ", "").length();
                                    if (spaces <= 5)
                                        addInfoMap.put(key, value.toString());
                                }
                            }
                            csb.setInfo(addInfoMap.keySet().stream().collect(Collectors.toMap(k -> k, k -> (List<String>) Arrays.asList(addInfoMap.get(k).toString()), (u,v) -> { throw new IllegalStateException(String.format("Duplicate key %s", u)); }, LinkedHashMap::new)));
                        }
            callSet = csb.build();
            listCallSet.add(callSet);
        }
        return SearchCallSetsResponse.newBuilder().setCallSets(listCallSet).setNextPageToken(nextPageToken).build();
    }

    @Override
    public SearchReferenceSetsResponse searchReferenceSets(SearchReferenceSetsRequest srsr) throws AvroRemoteException {
        List<String> list = new ArrayList<>();

        List<String> listModules = new ArrayList<>(MongoTemplateManager.getAvailableModules());
        Collections.sort(listModules);
        int start;
        int end;
        int pageSize;
        int pageToken = 0;
        String nextPageToken;

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

        List<ReferenceSet> listRef = Collections.synchronizedList(new ArrayList<>());
        List<String> modulesToReturn = start != 0 || end != listModules.size() ? listModules.subList(start, end) : listModules;
	    if (!modulesToReturn.isEmpty()) {
	        List<Collection<String>> splitModuleCollections = Helper.evenlySplitCollection(modulesToReturn, Runtime.getRuntime().availableProcessors() * 2);
	        ExecutorService executor = Executors.newFixedThreadPool(splitModuleCollections.size());
	        for (int i=0; i<splitModuleCollections.size(); i++) {
	            Collection<String> modules = splitModuleCollections.get(i);
	            Thread t = new Thread() {
	                public void run() {
	                    // add a Reference Set for each existing module
	                    for (String module : modules) {
	                        MongoTemplate mongoTemplate = MongoTemplateManager.get(module);
	                        String taxon = MongoTemplateManager.getTaxonName(module);
	                        String species = MongoTemplateManager.getSpecies(module);
	                        String taxoDesc = (species != null ? "Species: " + species : "") + (taxon != null && !taxon.equals(species) ? (species != null ? " ; " : "") + "Taxon: " + taxon : "");
	                        ReferenceSet referenceSet = ReferenceSet.newBuilder()
	                            .setId(module)
	                            .setName(module)
	                            .setMd5checksum("")    /* not supported at the time */
	                            .setSourceAccessions(list)
	                            .setNcbiTaxonId(MongoTemplateManager.getTaxonId(module))
	                            .setDescription(    (taxoDesc.isEmpty() ? "" : (taxoDesc + " ; "))
	                                                + mongoTemplate.getCollection(mongoTemplate.getCollectionName(GenotypingProject.class)).distinct(GenotypingProject.FIELDNAME_SEQUENCES, String.class).into(new ArrayList<>()).size() + " references ; "
	                                                + mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class)).estimatedDocumentCount() + " markers")
	                            .build();
	                        listRef.add(referenceSet);
	                    }
	                }
	            };
	            executor.execute(t);
	        }
	        executor.shutdown();
	        try {
	            executor.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);
	        } catch (InterruptedException e) {
	            throw new AvroRemoteException(e);
	        }
	    }

        return SearchReferenceSetsResponse.newBuilder().setReferenceSets(listRef).setNextPageToken(nextPageToken).build();
    }

    @Override
    public SearchVariantSetsResponse searchVariantSets(SearchVariantSetsRequest svsr) throws AvroRemoteException {

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
                vsmd.setKey(Constants.DESCRIPTION);
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
    public GigwaSearchVariantsResponse searchVariants(SearchVariantsRequest svr) throws AvroRemoteException {

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

            MongoCursor<Document> cursor = null;
            try
            {
                if (doSearch) {
                    // create a temp collection to store the result of the request
                    count = findVariants(gsvr);
                }
                else if (doCount || doBrowse) {
                    count = countVariants(gsvr, doBrowse);
                }

                if (count > 0 && doBrowse) {
                    String token = tokenManager.readToken(gsvr.getRequest());

                    if (token == null)
                        return null;

                    MongoTemplate mongoTemplate = MongoTemplateManager.get(info[0]);

                    MongoCollection<Document> tempVarColl = getTemporaryVariantCollection(info[0], token, false);
                    FindIterable<Document> iterable;
                    if (gsvr.getSelectedVariantIds() != null && tempVarColl.countDocuments() > 0) {
                        iterable = tempVarColl.find(); //when searching on variant IDs, retrieving all temporary collection
                    } else {
                        Collection<BasicDBList> variantQueryDBListCol = buildVariantDataQuery(gsvr, getSequenceIDsBeingFilteredOn(gsvr.getRequest().getSession(), info[0]), true);
                        //in this case, there is only one variantQueryDBList (no filtering on variant ids)
                        BasicDBList variantQueryDBList = !variantQueryDBListCol.isEmpty() ? variantQueryDBListCol.iterator().next() : new BasicDBList();

                        MongoCollection<Document> varCollForBuildingRows = tempVarColl.countDocuments() == 0 ? mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantData.class)) : tempVarColl;
                        iterable = varCollForBuildingRows.find(!variantQueryDBList.isEmpty() ? new BasicDBObject("$and", variantQueryDBList) : new BasicDBObject());
                    }
                    if(!parseBoolean(documentDbCompatMode))
                        iterable.collation(IExportHandler.collationObj);

                    if (gsvr.getSortBy() != null && gsvr.getSortBy().length() > 0)
                        iterable.sort(new BasicDBObject(gsvr.getSortBy(), Integer.valueOf("DESC".equalsIgnoreCase(gsvr.getSortDir()) ? -1 : 1)));
                    else if (mongoTemplate.findOne(new Query(new Criteria().andOperator(Criteria.where("_id").is(projId), Criteria.where(GenotypingProject.FIELDNAME_SEQUENCES + ".0").exists(true))), GenotypingProject.class) != null)
                        iterable.sort(new Document(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, 1).append(AbstractVariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, 1));
                    else
                        iterable.sort(new Document("_id", 1));  // no positions available in this project: let's sort variants by ID

                    iterable.skip(Integer.parseInt(gsvr.getPageToken()) * gsvr.getPageSize()).limit(gsvr.getPageSize());    // skip the results we don't want

                    cursor = iterable.iterator();
                }
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
                Collection<GenotypingSample> samples = new ArrayList<>();
                if (getGT) {
                    try {
                        samples = MgdbDao.getSamplesForProject(module, projId, gsvr.getCallSetIds().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(GigwaGa4ghServiceImpl.ID_SEPARATOR))).collect(Collectors.toList()));
                    } catch (ObjectNotFoundException ex) {
                        java.util.logging.Logger.getLogger(GigwaGa4ghServiceImpl.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }

                List<Variant> listVar = getVariantListFromDBCursor(module, Integer.parseInt(info[1]), cursor, samples);
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
    public SearchReferencesResponse searchReferences(SearchReferencesRequest srr) throws AvroRemoteException {

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
            List<String> accessions = new ArrayList<>();
            Map<String, Integer> mapSeq = new TreeMap<>(new AlphaNumericComparator());

            // allow search on checksum
            // but here only on one checksum v0.6.1 ?
            if (gsr.getMd5checksum() != null) {
                BasicDBObject query = new BasicDBObject();
                query.put(Sequence.FIELDNAME_CHECKSUM, gsr.getMd5checksum());
//                Document seq = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).find(query).first();
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
            ArrayList<BasicDBObject> pipeline = new ArrayList<>();
            pipeline.add(new BasicDBObject("$match", new BasicDBObject("_id", new BasicDBObject("$in", mapSeq.keySet()))));
            MongoCursor<Document> sqCursor = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(Sequence.class)).aggregate(pipeline).allowDiskUse(isAggregationAllowedToUseDisk()).iterator();

            Iterator<String> iteratorName = mapSeq.keySet().iterator();
            Iterator<Integer> iteratorId = mapSeq.values().iterator();

            // create and add the corresponding Reference for each sequence
            for (int i = start; i < end; i++) {
                Document sequence = !sqCursor.hasNext() ? null : sqCursor.next();
                String name = iteratorName.next();
                String projectId = Integer.toString(iteratorId.next());

                // Checksum : MD5 of the upper-case sequence excluding all whitespace characters (we usually don't have it)
                String md5 = sequence == null ? Helper.convertToMD5("") : (String) sequence.get(Sequence.FIELDNAME_CHECKSUM);
                Reference reference = Reference.newBuilder().setId(createId(module, projectId, name))
                        .setMd5checksum(md5 == null ? Helper.convertToMD5("") : md5)
                        .setName(name)
                        .setLength(sequence == null ? 0 : (long) sequence.get(Sequence.FIELDNAME_LENGTH))    // length == 0 when we don't have this information
                        .setSourceAccessions(accessions)
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
            varAnnField.put(VariantData.FIELDNAME_KNOWN_ALLELES, 1);
            varAnnField.put(VariantData.SECTION_ADDITIONAL_INFO, 1);
            Document variantRunDataObj = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(VariantRunData.class)).find(queryVarAnn).projection(varAnnField)
                .sort(new BasicDBObject(AbstractVariantData.SECTION_ADDITIONAL_INFO + "." + VariantRunData.FIELDNAME_ADDITIONAL_INFO_EFFECT_NAME, -1))  /*FIXME: this method should be called separately for each run*/
                .first();
            Document variantAnnotationObj = variantRunDataObj != null ? (Document) variantRunDataObj.get(VariantRunData.SECTION_ADDITIONAL_INFO) : null;
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

                    MongoCollection<Document> vcfHeaderColl = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class));
                    BasicDBList vcfHeaderQueryOrList = new BasicDBList();
                    for (String key : fieldHeader.keySet())
                        vcfHeaderQueryOrList.add(new BasicDBObject(key, new BasicDBObject("$exists", true)));

                    Document vcfHeaderEff = vcfHeaderColl.find(new BasicDBObject("$or", vcfHeaderQueryOrList)).projection(fieldHeader).first();

                    ArrayList<String> headerList = new ArrayList<>();
                    LinkedHashSet<String> usedHeaderSet = new LinkedHashSet<>();
                    if (!fAnnStyle)
                        headerList.add("Consequence");    // EFF style annotations
                    Document annInfo = (Document) ((Document) vcfHeaderEff.get(Constants.INFO_META_DATA)).get(fAnnStyle ? VcfImport.ANNOTATION_FIELDNAME_ANN : VcfImport.ANNOTATION_FIELDNAME_EFF);
                    if (annInfo == null && fAnnStyle)
                        annInfo = (Document) ((Document) vcfHeaderEff.get(Constants.INFO_META_DATA)).get(VcfImport.ANNOTATION_FIELDNAME_CSQ);
                    if (annInfo != null) {
                        header = (String) annInfo.get(Constants.DESCRIPTION);
                        if (header != null) {
                            // consider using the headers for additional info keySet
                            String sBeforeFieldList = fAnnStyle ? ": " : " (";
                            headerField = header.substring(header.indexOf(sBeforeFieldList) + sBeforeFieldList.length(), fAnnStyle ? header.length() : header.indexOf(")")).replaceAll("'", "").split("\\|");
                            for (String head : headerField)
                                headerList.add(head.replace("[", "").replace("]", "").trim());
                        }
                    }

                    List<AnalysisResult> listAnalysisResults = new ArrayList<>();

                    for (int i=0; i<tableTranscriptEffect.length; i++) {
                        ArrayList<String> values = new ArrayList<>();

                        if (!fAnnStyle) {    // EFF style annotations
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
                                            splitVals = value.split("-");    // sometimes used as separator
                                        try
                                        {
                                            allLocBuilder.setStart(Integer.parseInt(splitVals[0]));
                                        }
                                        catch (NumberFormatException ignored)
                                        {}

                                        if (allLocBuilder.getStart() == 0)
                                            continue;

                                        boolean fWorkingOnProtein = "Protein_position".equals(positionHeader);

                                        String sRefAllele = ((List<String>) variantRunDataObj.get(VariantData.FIELDNAME_KNOWN_ALLELES)).get(0);
                                        if (!fWorkingOnProtein)
                                            allLocBuilder.setEnd(allLocBuilder.getStart() + sRefAllele.length() - 1);
//                                        else
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
        Document result = MongoTemplateManager.get(module).getCollection(MongoTemplateManager.getMongoCollectionName(DBVCFHeader.class)).find(queryVarAnn).projection(varAnnField).first();
        if (result != null) {
            Document metaDataHeader = (Document) result.get(Constants.INFO_META_DATA);
            for (String key : metaDataHeader.keySet()) {
                annHeaders.put(key, (String) ((Document) metaDataHeader.get(key)).get(Constants.DESCRIPTION));
            }
            Document formatHeader = (Document) result.get(Constants.INFO_FORMAT_META_DATA);
            for (String key : formatHeader.keySet()) {
                annHeaders.put(key, (String) ((Document) formatHeader.get(key)).get(Constants.DESCRIPTION));
            }
        }
        return annHeaders;
    }

    @Override
    public TreeMap<String, HashMap<String, String>> getExportFormats() {
        TreeMap<String, HashMap<String, String>> exportFormats = new TreeMap<>();
        try {
            for (IExportHandler exportHandler : Stream.of(AbstractIndividualOrientedExportHandler.getIndividualOrientedExportHandlers().values(), AbstractMarkerOrientedExportHandler.getMarkerOrientedExportHandlers().values()).flatMap(Collection::stream).collect(Collectors.toList())) {
                HashMap<String, String> info = new HashMap<>();
                info.put("desc", exportHandler.getExportFormatDescription());
                info.put("supportedPloidyLevels", StringUtils.join(ArrayUtils.toObject(exportHandler.getSupportedPloidyLevels()), ";"));
                info.put("dataFileExtensions", StringUtils.join(exportHandler.getExportDataFileExtensions(), ";"));
                info.put("supportedVariantTypes", StringUtils.join(exportHandler.getSupportedVariantTypes(), ";"));
                exportFormats.put(exportHandler.getExportFormatName(), info);
            }
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException | SecurityException ex) {
            LOG.debug("error", ex);
        }
        return exportFormats;
    }

    public List<Comparable> searchVariantsLookup(String module, int projectId, String lookupText) throws AvroRemoteException {

        MongoTemplate mongoTemplate = MongoTemplateManager.get(module);

        List<Comparable> values = new ArrayList<>();

        MongoCollection<Document> collection = mongoTemplate.getCollection(MongoTemplateManager.getMongoCollectionName(VariantData.class));

        BasicDBObject whereQuery = new BasicDBObject();
        whereQuery.put("_id", Pattern.compile(".*\\Q" + lookupText + "\\E.*", Pattern.CASE_INSENSITIVE));

        int maxSize = 50;
        try {
            String variantIdLookupMaxSize = appConfig.get("variantIdLookupMaxSize");
            maxSize = Integer.parseInt(variantIdLookupMaxSize);
        } catch (Exception e) {
            LOG.debug("can't read variantIdLookupMaxSize in config, using maxSize=50");
        }

        MongoCursor<Document> cursor = collection.aggregate(
            Arrays.asList(
                Aggregates.match(whereQuery),
                Aggregates.group("$_id"),
                Aggregates.limit(maxSize+1)
            )
        ).iterator();

        try {
            while (cursor.hasNext()) {
                values.add((Comparable) cursor.next().get("_id"));
            }
        } finally {
           cursor.close();
        }

        if (values.size() > maxSize)
            return Arrays.asList("Too many results, please refine search!");

        return values;
    }
}