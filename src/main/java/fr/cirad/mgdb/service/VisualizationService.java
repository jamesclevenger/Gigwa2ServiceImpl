package fr.cirad.mgdb.service;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.log4j.Logger;
import org.bson.Document;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.data.mongodb.core.MongoTemplate;
import org.springframework.data.mongodb.core.query.Criteria;
import org.springframework.data.mongodb.core.query.Query;
import org.springframework.stereotype.Component;

import com.mongodb.BasicDBList;
import com.mongodb.BasicDBObject;
import com.mongodb.client.AggregateIterable;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;

import fr.cirad.mgdb.exporting.IExportHandler;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingProject;
import fr.cirad.mgdb.model.mongo.maintypes.GenotypingSample;
import fr.cirad.mgdb.model.mongo.maintypes.VariantData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData;
import fr.cirad.mgdb.model.mongo.maintypes.VariantRunData.VariantRunDataId;
import fr.cirad.mgdb.model.mongo.subtypes.ReferencePosition;
import fr.cirad.mgdb.model.mongo.subtypes.SampleGenotype;
import fr.cirad.mgdb.model.mongodao.MgdbDao;
import fr.cirad.model.GigwaDensityRequest;
import fr.cirad.model.GigwaSearchVariantsRequest;
import fr.cirad.model.GigwaVcfFieldPlotRequest;
import fr.cirad.tools.Helper;
import fr.cirad.tools.ProgressIndicator;
import fr.cirad.tools.mgdb.GenotypingDataQueryBuilder;
import fr.cirad.tools.mongo.MongoTemplateManager;
import fr.cirad.tools.security.base.AbstractTokenManager;
import fr.cirad.utils.Constants;

/**
 * A service class responsible for generating chart data
 * 
 * @author sempere
 */

@Component
public class VisualizationService {
    protected static final Logger LOG = Logger.getLogger(VisualizationService.class);
    
    @Autowired private AbstractTokenManager tokenManager;
    
	@Autowired private GigwaGa4ghServiceImpl ga4ghService;
	
    public boolean findDefaultRangeMinMax(GigwaDensityRequest gsvdr, String collectionName, ProgressIndicator progress)
    {
        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gsvdr.getVariantSetId(), 2);
        String sModule = info[0];

		final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);

		BasicDBList matchAndList = new BasicDBList();
		matchAndList.add(new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, gsvdr.getDisplayedSequence()));
        if ((gsvdr.getStart() != null && gsvdr.getStart() != -1) || (gsvdr.getEnd() != null && gsvdr.getEnd() != -1)) {
            BasicDBObject posCrit = new BasicDBObject();
            if (gsvdr.getStart() != null && gsvdr.getStart() != -1)
                posCrit.put("$gte", gsvdr.getStart());
            if (gsvdr.getEnd() != null && gsvdr.getEnd() != -1)
                posCrit.put("$lte", gsvdr.getEnd());
            matchAndList.add(new BasicDBObject(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE, posCrit));
        }
		if (gsvdr.getDisplayedVariantType() != null)
			matchAndList.add(new BasicDBObject(VariantData.FIELDNAME_TYPE, gsvdr.getDisplayedVariantType()));
		BasicDBObject match = new BasicDBObject("$match", new BasicDBObject("$and", matchAndList));

		String startFieldPath = VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE;
		BasicDBObject sort = new BasicDBObject("$sort", new BasicDBObject(startFieldPath, 1));
		BasicDBObject limit = new BasicDBObject("$limit", 1);
		MongoCursor<Document> cursor = mongoTemplate.getCollection(collectionName).aggregate(Arrays.asList(match, sort, limit)).iterator();
		if (!cursor.hasNext()) {
			if (progress != null)
				progress.markAsComplete();
			return false;	// no variant found matching filter
		}
		Document aggResult = (Document) cursor.next();
		if (gsvdr.getDisplayedRangeMin() == null)
			gsvdr.setDisplayedRangeMin((Long) Helper.readPossiblyNestedField(aggResult, startFieldPath, "; "));

		sort = new BasicDBObject("$sort", new BasicDBObject(startFieldPath, -1));
		cursor = mongoTemplate.getCollection(collectionName).aggregate(Arrays.asList(match, sort, limit)).collation(IExportHandler.collationObj).iterator();
		if (!cursor.hasNext()) {
			if (progress != null)
				progress.markAsComplete();
			return false;	// no variant found matching filter
		}
		aggResult = (Document) cursor.next();
		if (gsvdr.getDisplayedRangeMax() == null)
			gsvdr.setDisplayedRangeMax((Long) Helper.readPossiblyNestedField(aggResult, startFieldPath, "; "));
		return true;
	}
    
    public Map<Long, Long> selectionDensity(GigwaDensityRequest gdr) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gdr.getVariantSetId(), 2);
        String sModule = info[0];

        ProgressIndicator progress = new ProgressIndicator(tokenManager.readToken(gdr.getRequest()), new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "variant density on sequence " + gdr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Collection<BasicDBList> variantQueryDBListColl = ga4ghService.buildVariantDataQuery(gdr, ga4ghService.getSequenceIDsBeingFilteredOn(gdr.getRequest().getSession(), sModule), true);
        final BasicDBList variantQueryDBList = variantQueryDBListColl.size() == 1 ? variantQueryDBListColl.iterator().next() : new BasicDBList();

        MongoCollection<Document> tmpVarColl = ga4ghService.getTemporaryVariantCollection(sModule, tokenManager.readToken(gdr.getRequest()), false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0)
        {
            progress.setError(Constants.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        final String mainVarCollName = mongoTemplate.getCollectionName(VariantData.class), usedVarCollName = nTempVarCount == 0 ? mainVarCollName : tmpVarColl.getNamespace().getCollectionName();
        final ConcurrentHashMap<Long, Long> result = new ConcurrentHashMap<Long, Long>();

        if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gdr, usedVarCollName, progress))
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
                        long partialCount = mongoTemplate.getCollection(usedVarCollName).countDocuments(new BasicDBObject("$and", queryList));
                        nTotalTreatedVariantCount.addAndGet((int) partialCount);
                        result.put(rangeMin + (chunkIndex*intervalSize), partialCount);
                        finalProgress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
                    }
                }
            };

            if (chunkIndex%GigwaGa4ghServiceImpl.INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS  == (GigwaGa4ghServiceImpl.INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS-1))
                t.run();    // run synchronously
            else
            {
                threadsToWaitFor.add(t);
                t.start();    // run asynchronously for better speed
            }
        }

        if (progress.isAborted())
                return null;

        for (Thread ttwf : threadsToWaitFor)    // wait for all threads before moving to next phase
                ttwf.join();

        progress.setCurrentStepProgress(100);
        LOG.debug("selectionDensity treated " + nTotalTreatedVariantCount.get() + " variants in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

		return new TreeMap<Long, Long>(result);
	}

    // TODO: Refactor this?
    private void mergeVariantQueryDBList(BasicDBObject matchStage, BasicDBList variantQueryDBList) {
    	Iterator<Object> queryItems = variantQueryDBList.iterator();
		while (queryItems.hasNext()) {
			BasicDBObject queryItem = (BasicDBObject)queryItems.next();
			for (String key : queryItem.keySet()) {
				if (queryItem.get(key) instanceof BasicDBObject) {
					BasicDBObject queryItemElement = (BasicDBObject)queryItem.get(key);
					if (matchStage.containsKey(key)) {
						if (matchStage.get(key) instanceof BasicDBObject) {
							BasicDBObject matchStageElement = (BasicDBObject)matchStage.get(key);
							for (String elementKey : queryItemElement.keySet()) {
								if (matchStageElement.containsKey(elementKey)) {
									if (elementKey.equals("$lt") || elementKey.equals("$lte")) {
										matchStageElement.put(elementKey, Math.min(matchStageElement.getLong(elementKey), queryItemElement.getLong(elementKey)));
									} else if (elementKey.equals("$gt") || elementKey.equals("$gte")) {
										matchStageElement.put(elementKey, Math.max(matchStageElement.getLong(elementKey), queryItemElement.getLong(elementKey)));
									} else {
										matchStageElement.put(elementKey, queryItemElement.get(elementKey));
									}
								} else {
									matchStageElement.put(elementKey, queryItemElement.get(elementKey));
								}
							}
						} else {
							matchStage.put(key, queryItemElement);
						}
					} else {
						matchStage.put(key, queryItemElement);
					}
				} else {
					matchStage.put(key, queryItem.get(key));
				}
			}
		}
    }

    public Map<Long, Double> selectionFst(GigwaDensityRequest gdr) throws Exception {
    	long before = System.currentTimeMillis();

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gdr.getVariantSetId(), 2);
        String sModule = info[0];

		ProgressIndicator progress = new ProgressIndicator(tokenManager.readToken(gdr.getRequest()), new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "Fst estimate on sequence " + gdr.getDisplayedSequence()});
		ProgressIndicator.registerProgressIndicator(progress);

		final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        final BasicDBList variantQueryDBList = ga4ghService.buildVariantDataQuery(gdr, ga4ghService.getSequenceIDsBeingFilteredOn(gdr.getRequest().getSession(), sModule), true).iterator().next();   // there's only one in this case

		MongoCollection<Document> tmpVarColl = ga4ghService.getTemporaryVariantCollection(sModule, tokenManager.readToken(gdr.getRequest()), false);
		long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
		if (GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0)
		{
			progress.setError(Constants.MESSAGE_TEMP_RECORDS_NOT_FOUND);
			return null;
		}

		final String vrdCollName = mongoTemplate.getCollectionName(VariantRunData.class);
		final boolean useTempColl = (nTempVarCount != 0);
		final String usedVarCollName = useTempColl ? tmpVarColl.getNamespace().getCollectionName() : vrdCollName;
		final ConcurrentHashMap<Long, Double> result = new ConcurrentHashMap<Long, Double>();

		if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
			if (!findDefaultRangeMinMax(gdr, usedVarCollName, progress))
				return result;

		final AtomicInteger nTotalTreatedVariantCount = new AtomicInteger(0);
		final int intervalSize = Math.max(1, (int) ((gdr.getDisplayedRangeMax() - gdr.getDisplayedRangeMin()) / gdr.getDisplayedRangeIntervalCount()));
		final long rangeMin = gdr.getDisplayedRangeMin();
		final ProgressIndicator finalProgress = progress;

		List<BasicDBObject> baseQuery = buildFstQuery(gdr, useTempColl);

		int nConcurrentThreads = Math.min(Runtime.getRuntime().availableProcessors(), GigwaGa4ghServiceImpl.INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS);
		ExecutorService executor = Executors.newFixedThreadPool(nConcurrentThreads);

		for (int i=0; i<gdr.getDisplayedRangeIntervalCount(); i++) {
			BasicDBObject initialMatchStage = new BasicDBObject();
			initialMatchStage.put(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, gdr.getDisplayedSequence());
			if (gdr.getDisplayedVariantType() != null)
				initialMatchStage.put(VariantData.FIELDNAME_TYPE, gdr.getDisplayedVariantType());
			BasicDBObject positionSettings = new BasicDBObject();
			positionSettings.put("$gte", gdr.getDisplayedRangeMin() + (i*intervalSize));
			positionSettings.put(i < gdr.getDisplayedRangeIntervalCount() - 1 ? "$lt" : "$lte", i < gdr.getDisplayedRangeIntervalCount() - 1 ? gdr.getDisplayedRangeMin() + ((i+1)*intervalSize) : gdr.getDisplayedRangeMax());
			String startSitePath = VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE;
			initialMatchStage.put(startSitePath, positionSettings);
			if (nTempVarCount == 0 && !variantQueryDBList.isEmpty())
				mergeVariantQueryDBList(initialMatchStage, variantQueryDBList);
			final long chunkIndex = i;

			List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
			windowQuery.set(0, new BasicDBObject("$match", initialMatchStage));

			//try { System.out.println(new ObjectMapper().writeValueAsString(windowQuery)); }
            //catch (Exception ignored) {}

            Thread t = new Thread() {
            	public void run() {
            		if (finalProgress.isAborted())
            			return;

            		Iterator<Document> it = mongoTemplate.getCollection(usedVarCollName).aggregate(windowQuery).allowDiskUse(ga4ghService.isAggregationAllowedToUseDisk()).iterator();

            		if (finalProgress.isAborted())
            			return;

            		/* Structure of a resulting document : {
        			 * 		_id: ...,
        			 * 		alleleMax: 1,
        			 * 		populations: [
        			 * 			{sampleSize: 100, alleles: [
        			 * 				{allele: 0, alleleFrequency: 0.45, heterozygoteFrequency: 0.31},
        			 * 				{allele: 1, alleleFrequency: 0.55, heterozygoteFrequency: 0.31},
        			 * 			]},
        			 * 			{...}
        			 * 		]
        			 * }
        			 */

        			double weightedFstSum = 0;
        			double fstWeight = 0;

        			while (it.hasNext()) {
        				Document variantResult = it.next();
        				//String variantId = variantResult.getString("_id");

        				List<Document> populations = variantResult.getList(FST_RES_POPULATIONS, Document.class);
        				if (populations.size() < 2) {
        					// Can not compute Fst with a single population
        					// One of the populations has no valid data
        					continue;
        				}
        				int numPopulations = populations.size();  // r : Number of samples to consider
        				int numAlleles = variantResult.getInteger(FST_RES_ALLELEMAX) + 1;

        				// Transposition to [allele][sample] instead of the original [sample][allele] is important to simplify further computations
        				int[] sampleSizes = new int[numPopulations];  // n_i = sampleSizes[population] : Size of the population samples (with missing data filtered out)
        				double[][] alleleFrequencies = new double[numAlleles][numPopulations];  // p_i = alleleFrequencies[allele][population] : Allele frequency in the given population
        				double[][] hetFrequencies = new double[numAlleles][numPopulations];  // h_i = hetFrequencies[allele][population] : Proportion of heterozygotes with the given allele in the given population
        				//double[] averageAlleleFrequencies = new double[numAlleles];  // p¯ = averageAlleleFrequencies[allele] : Average frequency of the allele over all populations
    					//double[] alleleVariance = new double[numAlleles];  // s² = alleleVariance[allele] : Variance of the allele frequency over the populations
    					//double[] averageHetFrequencies = new double[numAlleles];  // h¯ = averageHetFrequencies : Proportion of heterozygotes with the given allele over all populations

    					//Arrays.fill(averageAlleleFrequencies, 0);
    					//Arrays.fill(alleleVariance, 0);
    					//Arrays.fill(averageHetFrequencies, 0);

    					for (int allele = 0; allele < numAlleles; allele++) {
    						Arrays.fill(alleleFrequencies[allele], 0);
    						Arrays.fill(hetFrequencies[allele], 0);
    					}

    					int popIndex = 0;
        				for (Document populationResult : populations) {
        					int sampleSize = populationResult.getInteger(FST_RES_SAMPLESIZE);
        					List<Document> alleles = populationResult.getList(FST_RES_ALLELES, Document.class);

        					for (Document alleleResult : alleles) {
        						int allele = alleleResult.getInteger(FST_RES_ALLELEID);
        						alleleFrequencies[allele][popIndex] = alleleResult.getDouble(FST_RES_ALLELEFREQUENCY);
        						hetFrequencies[allele][popIndex] = alleleResult.getDouble(FST_RES_HETEROZYGOTEFREQUENCY);
        					}

        					sampleSizes[popIndex] = sampleSize;
        					popIndex += 1;
        				}

        				double averageSampleSize = (double)IntStream.of(sampleSizes).sum() / numPopulations;  // n¯ : Average sample size
        				double totalSize = averageSampleSize * numPopulations;  // r × n¯
        				double sampleSizeCorrection = (totalSize - IntStream.of(sampleSizes).mapToDouble(size -> size*size / totalSize).sum() / (numPopulations - 1));  // n_c

        				for (int allele = 0; allele < numAlleles; allele++) {
        					// Compute weighted averages of allele frequencies (p¯) and heterozygote proportions (h¯)
        					double averageAlleleFrequency = 0.0;
        					double averageHetFrequency = 0.0;
        					for (popIndex = 0; popIndex < numPopulations; popIndex++) {
        						averageAlleleFrequency += sampleSizes[popIndex] * alleleFrequencies[allele][popIndex] / totalSize;
        						averageHetFrequency += sampleSizes[popIndex] * hetFrequencies[allele][popIndex] / totalSize;
        					}

        					// Compute allele frequency variance (s²)
        					double alleleVariance = 0.0;
        					for (popIndex = 0; popIndex < numPopulations; popIndex++) {
        						alleleVariance += sampleSizes[popIndex] * Math.pow(alleleFrequencies[allele][popIndex] - averageAlleleFrequency, 2) / (averageSampleSize * (numPopulations - 1));
        					}

        					// a = (n¯/nc) × (s² - (1 / (n¯-1))(p¯(1-p¯) - s²(r-1) / r - h¯/4))
							double populationVariance = (
									(averageSampleSize / sampleSizeCorrection) * (
										alleleVariance - (1 / (averageSampleSize - 1)) * (
											averageAlleleFrequency * (1 - averageAlleleFrequency) -
											(numPopulations - 1) * alleleVariance / numPopulations -
											averageHetFrequency / 4
										)
									)
								);

							// b = (n¯ / (n¯-1))(p¯(1-p¯) - s²(r-1) / r - h¯(2n¯-1)/4n¯)
							double individualVariance = (averageSampleSize / (averageSampleSize - 1)) * (
									averageAlleleFrequency * (1 - averageAlleleFrequency) -
									alleleVariance * (numPopulations-1) / numPopulations -
									averageHetFrequency * (2*averageSampleSize - 1) / (4 * averageSampleSize)
								);

							// c = h¯/2
							double gameteVariance = averageHetFrequency / 2;

							if (!Double.isNaN(populationVariance) && !Double.isNaN(individualVariance) && !Double.isNaN(gameteVariance)) {
								weightedFstSum += populationVariance;
								fstWeight += populationVariance + individualVariance + gameteVariance;
							}
        				}
        			}

        			result.put(rangeMin + (chunkIndex*intervalSize), weightedFstSum / fstWeight);
        			finalProgress.setCurrentStepProgress((short) result.size() * 100 / gdr.getDisplayedRangeIntervalCount());
        		}
            };

           executor.execute(t);
		}

		executor.shutdown();
		executor.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);

		if (progress.isAborted())
			return null;

		progress.setCurrentStepProgress(100);
		LOG.info("selectionDensity treated " + nTotalTreatedVariantCount.get() + " variants in " + (System.currentTimeMillis() - before)/1000f + "s");
		progress.markAsComplete();

		return new TreeMap<Long, Double>(result);
    }
    
    public List<Map<Long, Double>> selectionTajimaD(GigwaDensityRequest gdr) throws Exception {
		long before = System.currentTimeMillis();

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gdr.getVariantSetId(), 2);
        String sModule = info[0];

		ProgressIndicator progress = new ProgressIndicator(tokenManager.readToken(gdr.getRequest()), new String[] {"Calculating " + (gdr.getDisplayedVariantType() != null ? gdr.getDisplayedVariantType() + " " : "") + "Tajima's D on sequence " + gdr.getDisplayedSequence()});
		ProgressIndicator.registerProgressIndicator(progress);

		final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
		final BasicDBList variantQueryDBList = ga4ghService.buildVariantDataQuery(gdr, ga4ghService.getSequenceIDsBeingFilteredOn(gdr.getRequest().getSession(), sModule), true).iterator().next();   // there's only one in this case

		MongoCollection<Document> tmpVarColl = ga4ghService.getTemporaryVariantCollection(sModule, tokenManager.readToken(gdr.getRequest()), false);
		long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
		if (GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gdr, false).size() > 0 && nTempVarCount == 0) {
			progress.setError(Constants.MESSAGE_TEMP_RECORDS_NOT_FOUND);
			return null;
		}

		final String vrdCollName = mongoTemplate.getCollectionName(VariantRunData.class);
		final boolean useTempColl = (nTempVarCount != 0);
		final String usedVarCollName = useTempColl ? tmpVarColl.getNamespace().getCollectionName() : vrdCollName;
		final ConcurrentHashMap<Long, Double> tajimaD = new ConcurrentHashMap<Long, Double>();
		final ConcurrentHashMap<Long, Double> segregatingSites = new ConcurrentHashMap<Long, Double>();
		final List<Map<Long, Double>> result = Arrays.asList(tajimaD, segregatingSites);

		if (gdr.getDisplayedRangeMin() == null || gdr.getDisplayedRangeMax() == null)
			if (!findDefaultRangeMinMax(gdr, usedVarCollName, progress))
				return result;

		List<BasicDBObject> baseQuery = buildTajimaDQuery(gdr, useTempColl);

		int nConcurrentThreads = Math.min(Runtime.getRuntime().availableProcessors(), GigwaGa4ghServiceImpl.INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS);
		final int intervalSize = Math.max(1, (int) ((gdr.getDisplayedRangeMax() - gdr.getDisplayedRangeMin()) / gdr.getDisplayedRangeIntervalCount()));
		ExecutorService executor = Executors.newFixedThreadPool(nConcurrentThreads);

		for (int i=0; i<gdr.getDisplayedRangeIntervalCount(); i++) {
			BasicDBObject initialMatchStage = new BasicDBObject();
			initialMatchStage.put(VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_SEQUENCE, gdr.getDisplayedSequence());
			if (gdr.getDisplayedVariantType() != null)
				initialMatchStage.put(VariantData.FIELDNAME_TYPE, gdr.getDisplayedVariantType());
			BasicDBObject positionSettings = new BasicDBObject();
			positionSettings.put("$gte", gdr.getDisplayedRangeMin() + (i*intervalSize));
			positionSettings.put(i < gdr.getDisplayedRangeIntervalCount() - 1 ? "$lt" : "$lte", i < gdr.getDisplayedRangeIntervalCount() - 1 ? gdr.getDisplayedRangeMin() + ((i+1)*intervalSize) : gdr.getDisplayedRangeMax());
			String startSitePath = VariantData.FIELDNAME_REFERENCE_POSITION + "." + ReferencePosition.FIELDNAME_START_SITE;
			initialMatchStage.put(startSitePath, positionSettings);
			if (nTempVarCount == 0 && !variantQueryDBList.isEmpty())
				mergeVariantQueryDBList(initialMatchStage, variantQueryDBList);
			final long chunkIndex = i;

			List<BasicDBObject> windowQuery = new ArrayList<BasicDBObject>(baseQuery);
			windowQuery.set(0, new BasicDBObject("$match", initialMatchStage));

			Thread t = new Thread() {
				
				public void run() {
					if (progress.isAborted())
						return;

					AggregateIterable<Document> queryResult = mongoTemplate.getCollection(usedVarCollName).aggregate(windowQuery).allowDiskUse(ga4ghService.isAggregationAllowedToUseDisk());

					if (progress.isAborted())
						return;

					Document chunk = queryResult.first();  // There's only one interval per query

					long intervalStart = gdr.getDisplayedRangeMin() + (chunkIndex*intervalSize);
					if (chunk != null) {
						double value = chunk.getDouble(TJD_RES_TAJIMAD);
						double sites = (double)chunk.getInteger(TJD_RES_SEGREGATINGSITES);
						segregatingSites.put(intervalStart, sites);
						tajimaD.put(intervalStart, value);
					} else {
						segregatingSites.put(intervalStart, 0.0);
						tajimaD.put(intervalStart, Double.NaN);
					}
        			progress.setCurrentStepProgress((short) segregatingSites.size() * 100 / gdr.getDisplayedRangeIntervalCount());
				}
			};

			executor.execute(t);
		}

		executor.shutdown();
		executor.awaitTermination(Integer.MAX_VALUE, TimeUnit.DAYS);

		if (progress.isAborted())
			return null;

		progress.setCurrentStepProgress(100);
		LOG.debug("selectionTajimaD treated in " + (System.currentTimeMillis() - before)/1000f + "s");
		progress.markAsComplete();

		return result;
	}

    private static final String GENOTYPE_DATA_S2_DATA = "dt";
    private static final String GENOTYPE_DATA_S5_SPKEYVAL = "sk";
    private static final String GENOTYPE_DATA_S8_SAMPLE = "sa";
    private static final String GENOTYPE_DATA_S7_SAMPLEID = "si";
    private static final String GENOTYPE_DATA_S7_GENOTYPE = "gy";
    private static final String GENOTYPE_DATA_S7_POSITION = "ss";
    private static final String GENOTYPE_DATA_S10_VARIANTID = "vi";
    private static final String GENOTYPE_DATA_S10_INDIVIDUALID = "ii";
    private static final String GENOTYPE_DATA_S10_SAMPLEINDEX = "sx";

    private List<BasicDBObject> buildGenotypeDataQuery(GigwaDensityRequest gdr, boolean useTempColl, Map<String, List<GenotypingSample>> individualToSampleListMap, boolean keepPosition) {
    	String info[] = GigwaSearchVariantsRequest.getInfoFromId(gdr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        GenotypingProject genotypingProject = mongoTemplate.findById(Integer.valueOf(projId), GenotypingProject.class);

        boolean fIsMultiRunProject = genotypingProject.getRuns().size() > 1;
        boolean fGotMultiSampleIndividuals = false;


        for (List<GenotypingSample> samplesForAGivenIndividual : individualToSampleListMap.values()) {
            if (samplesForAGivenIndividual.size() > 1) {
                fGotMultiSampleIndividuals = true;
                break;
            }
        }

    	List<BasicDBObject> pipeline = new ArrayList<BasicDBObject>();

    	// Stage 1 : placeholder for initial match stage
    	pipeline.add(null);

    	if (useTempColl) {
	    	// Stage 2 : Lookup from temp collection to variantRunData
	    	BasicDBObject lookup = new BasicDBObject();
	    	lookup.put("from", "variantRunData");
	    	lookup.put("localField", "_id");
	    	lookup.put("foreignField", "_id.vi");
	    	lookup.put("as", GENOTYPE_DATA_S2_DATA);
	    	pipeline.add(new BasicDBObject("$lookup", lookup));

	    	// Stage 3 : Unwind data
	    	pipeline.add(new BasicDBObject("$unwind", "$" + GENOTYPE_DATA_S2_DATA));

	    	// Stage 4 : Keep only the right project
	    	pipeline.add(new BasicDBObject("$match", new BasicDBObject(GENOTYPE_DATA_S2_DATA + "._id.pi", projId)));
    	} else {
    		// Stage 4 : Keep only the right project
    		pipeline.add(new BasicDBObject("$match", new BasicDBObject("_id.pi", projId)));
    	}

    	if (fGotMultiSampleIndividuals) {
    		if (useTempColl) {
	    		// Stage 5 : Convert samples to an array
	    		pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(GENOTYPE_DATA_S5_SPKEYVAL, new BasicDBObject("$objectToArray", "$" + GENOTYPE_DATA_S2_DATA + "." + VariantRunData.FIELDNAME_SAMPLEGENOTYPES))));
    		} else {
    			// Stage 5 : Convert samples to an array
	    		pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(GENOTYPE_DATA_S5_SPKEYVAL, new BasicDBObject("$objectToArray", "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES))));
    		}

    		// Stage 6 : Unwind samples
    		pipeline.add(new BasicDBObject("$unwind", "$" + GENOTYPE_DATA_S5_SPKEYVAL));

    		// Stage 7 : Convert key and value
    		BasicDBObject spkeyval = new BasicDBObject();
    		spkeyval.put(GENOTYPE_DATA_S7_SAMPLEID, new BasicDBObject("$toInt", "$" + GENOTYPE_DATA_S5_SPKEYVAL + ".k"));
    		spkeyval.put(GENOTYPE_DATA_S7_GENOTYPE, "$" + GENOTYPE_DATA_S5_SPKEYVAL + ".v.gt");
    		if (keepPosition)
    			spkeyval.put(GENOTYPE_DATA_S7_POSITION, "$" + VariantData.FIELDNAME_REFERENCE_POSITION + ".ss");
    		pipeline.add(new BasicDBObject("$project", spkeyval));

    		// Stage 8 : Lookup samples
    		BasicDBObject sampleLookup = new BasicDBObject();
    		sampleLookup.put("from", "samples");
    		sampleLookup.put("localField", GENOTYPE_DATA_S7_SAMPLEID);
    		sampleLookup.put("foreignField", "_id");
    		sampleLookup.put("as", GENOTYPE_DATA_S8_SAMPLE);
    		pipeline.add(new BasicDBObject("$lookup", sampleLookup));

    		// Stage 9 : Get first sample (shouldn't get more than one anyway)
    		BasicDBObject firstSample = new BasicDBObject();
    		firstSample.put(GENOTYPE_DATA_S7_GENOTYPE, 1);
    		if (keepPosition)
    			firstSample.put(GENOTYPE_DATA_S7_POSITION, 1);
    		firstSample.put(GENOTYPE_DATA_S8_SAMPLE, new BasicDBObject("$arrayElemAt", Arrays.asList("$" + GENOTYPE_DATA_S8_SAMPLE, 0)));
    		pipeline.add(new BasicDBObject("$project", firstSample));

    		// Stage 10 : Regroup individual runs
    		BasicDBObject individualGroup = new BasicDBObject();
    		BasicDBObject individualGroupId = new BasicDBObject();
    		individualGroupId.put(GENOTYPE_DATA_S10_VARIANTID, "$_id");
    		individualGroupId.put(GENOTYPE_DATA_S10_INDIVIDUALID, "$" + GENOTYPE_DATA_S8_SAMPLE + "." + GenotypingSample.FIELDNAME_INDIVIDUAL);
    		individualGroup.put("_id", individualGroupId);
    		individualGroup.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$addToSet", "$" + GENOTYPE_DATA_S7_GENOTYPE));
    		individualGroup.put(GENOTYPE_DATA_S10_SAMPLEINDEX, new BasicDBObject("$min", "$" + GENOTYPE_DATA_S8_SAMPLE + "._id"));
    		if (keepPosition)
    			individualGroup.put(GENOTYPE_DATA_S7_POSITION, new BasicDBObject("$first", "$" + GENOTYPE_DATA_S7_POSITION));
    		pipeline.add(new BasicDBObject("$group", individualGroup));

    		// Stage 11 : Weed out incoherent genotypes
    		pipeline.add(new BasicDBObject("$match", new BasicDBObject(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$size", 1))));

    		// Stage 12 : Group back by variant
    		BasicDBObject variantGroup = new BasicDBObject();
    		variantGroup.put("_id", "$_id." + GENOTYPE_DATA_S10_VARIANTID);
    		BasicDBObject spObject = new BasicDBObject();
    		spObject.put("k", new BasicDBObject("$toString", "$" + GENOTYPE_DATA_S10_SAMPLEINDEX));
    		spObject.put("v", new BasicDBObject(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$arrayElemAt", Arrays.asList("$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES, 0))));
    		variantGroup.put("sp", new BasicDBObject("$push", spObject));
    		if (keepPosition)
    			variantGroup.put(GENOTYPE_DATA_S7_POSITION, new BasicDBObject("$first", "$" + GENOTYPE_DATA_S7_POSITION));
    		pipeline.add(new BasicDBObject("$group", variantGroup));

    		// Stage 13 : Convert back to sp object
    		pipeline.add(new BasicDBObject("$project", new BasicDBObject(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$arrayToObject", "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES))));
    	} else {
    		if (useTempColl) {
		    	// Stage 5 : Group runs
		    	BasicDBObject groupRuns = new BasicDBObject();
		    	groupRuns.put("_id", "$_id");
		    	groupRuns.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$first", "$" + GENOTYPE_DATA_S2_DATA + "." + VariantRunData.FIELDNAME_SAMPLEGENOTYPES));
		    	pipeline.add(new BasicDBObject("$group", groupRuns));
    		} else if (fIsMultiRunProject) {
    			// Stage 5 : Group runs
		    	BasicDBObject groupRuns = new BasicDBObject();
		    	groupRuns.put("_id", "$_id");
		    	groupRuns.put(VariantRunData.FIELDNAME_SAMPLEGENOTYPES, new BasicDBObject("$first", "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES));
		    	pipeline.add(new BasicDBObject("$group", groupRuns));
    		}
    	}

    	return pipeline;
    }


    private static final String FST_S14_POPULATIONGENOTYPES = "pg";
    private static final String FST_S15_POPULATION = "pp";
    private static final String FST_S19_GENOTYPE = "gn";
    private static final String FST_S20_HETEROZYGOTE = "ht";
    private static final String FST_S22_VARIANTID = "vi";
    private static final String FST_S22_POPULATIONID = "pi";
    private static final String FST_S22_ALLELECOUNT = "ac";
    private static final String FST_S22_HETEROZYGOTECOUNT = "hc";

    private static final String FST_RES_SAMPLESIZE = "ss";
    private static final String FST_RES_ALLELEID = "al";
    private static final String FST_RES_ALLELEMAX = "am";
    private static final String FST_RES_ALLELEFREQUENCY = "af";
    private static final String FST_RES_HETEROZYGOTEFREQUENCY = "hf";
    private static final String FST_RES_ALLELES = "as";
    private static final String FST_RES_POPULATIONS = "ps";

    private List<BasicDBObject> buildFstQuery(GigwaDensityRequest gdr, boolean useTempColl) {
    	String info[] = GigwaSearchVariantsRequest.getInfoFromId(gdr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

    	List<Collection<String>> selectedIndividuals = new ArrayList<Collection<String>>();
        if (gdr.getDisplayedAdditionalGroups() == null) {
        	selectedIndividuals.add(gdr.getCallSetIds().size() == 0 ? MgdbDao.getProjectIndividuals(sModule, projId) : gdr.getCallSetIds().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(IGigwaService.ID_SEPARATOR))).collect(Collectors.toSet()));
        	selectedIndividuals.add(gdr.getCallSetIds2().size() == 0 ? MgdbDao.getProjectIndividuals(sModule, projId) : gdr.getCallSetIds2().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(IGigwaService.ID_SEPARATOR))).collect(Collectors.toSet()));
        } else {
        	for (Collection<String> group : gdr.getDisplayedAdditionalGroups()) {
        		selectedIndividuals.add(group.size() == 0 ? MgdbDao.getProjectIndividuals(sModule, projId) : group.stream().map(csi -> csi.substring(1 + csi.lastIndexOf(IGigwaService.ID_SEPARATOR))).collect(Collectors.toSet()));
        	}
        }
        TreeMap<String, List<GenotypingSample>> individualToSampleListMap = new TreeMap<String, List<GenotypingSample>>();
        for (Collection<String> group : selectedIndividuals) {
        	individualToSampleListMap.putAll(MgdbDao.getSamplesByIndividualForProject(sModule, projId, group));
        }

    	List<BasicDBObject> pipeline = buildGenotypeDataQuery(gdr, useTempColl, individualToSampleListMap, false);

    	// Stage 14 : Get populations genotypes
    	BasicDBList populationGenotypes = new BasicDBList();
    	for (Collection<String> group : selectedIndividuals) {
    		populationGenotypes.add(getFullPathToGenotypes(sModule, projId, group, individualToSampleListMap));
    	}

    	BasicDBObject projectGenotypes = new BasicDBObject(FST_S14_POPULATIONGENOTYPES, populationGenotypes);
    	pipeline.add(new BasicDBObject("$project", projectGenotypes));

    	// Stage 15 : Split by population
    	BasicDBObject unwindPopulations = new BasicDBObject();
    	unwindPopulations.put("path", "$" + FST_S14_POPULATIONGENOTYPES);
    	unwindPopulations.put("includeArrayIndex", FST_S15_POPULATION);
    	pipeline.add(new BasicDBObject("$unwind", unwindPopulations));

    	// Stage 16 : Compute sample size
    	BasicDBObject sampleSizeMapping = new BasicDBObject();
    	sampleSizeMapping.put("input", "$" + FST_S14_POPULATIONGENOTYPES);
    	sampleSizeMapping.put("in", new BasicDBObject("$cmp", Arrays.asList("$$this", null)));
    	BasicDBObject addSampleSize = new BasicDBObject("$sum", new BasicDBObject("$map", sampleSizeMapping));
    	pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(FST_RES_SAMPLESIZE, addSampleSize)));

    	// Stage 17 : Unwind by genotype
    	pipeline.add(new BasicDBObject("$unwind", "$" + FST_S14_POPULATIONGENOTYPES));

    	// Stage 18 : Eliminate missing genotypes
    	BasicDBObject matchMissing = new BasicDBObject(FST_S14_POPULATIONGENOTYPES, new BasicDBObject("$ne", null));
    	pipeline.add(new BasicDBObject("$match", matchMissing));

    	// Stage 19 : Split genotype strings
    	BasicDBObject projectSplitGenotypes = new BasicDBObject();
    	projectSplitGenotypes.put(FST_S15_POPULATION, 1);
    	projectSplitGenotypes.put(FST_RES_SAMPLESIZE, 1);
    	BasicDBObject splitMapping = new BasicDBObject();
    	splitMapping.put("input", new BasicDBObject("$split", Arrays.asList("$" + FST_S14_POPULATIONGENOTYPES, "/")));
    	splitMapping.put("in", new BasicDBObject("$toInt", "$$this"));
    	projectSplitGenotypes.put(FST_S19_GENOTYPE, new BasicDBObject("$map", splitMapping));
    	pipeline.add(new BasicDBObject("$project", projectSplitGenotypes));

    	// Stage 20 : Detect heterozygotes
    	BasicDBList genotypeElements = new BasicDBList();  // TODO : Ploidy ?
    	genotypeElements.add(new BasicDBObject("$arrayElemAt", Arrays.asList("$" + FST_S19_GENOTYPE, 0)));
    	genotypeElements.add(new BasicDBObject("$arrayElemAt", Arrays.asList("$" + FST_S19_GENOTYPE, 1)));
    	BasicDBObject addHeterozygote = new BasicDBObject(FST_S20_HETEROZYGOTE, new BasicDBObject("$ne", genotypeElements));
    	pipeline.add(new BasicDBObject("$addFields", addHeterozygote));

    	// Stage 21 : Unwind alleles
    	pipeline.add(new BasicDBObject("$unwind", "$" + FST_S19_GENOTYPE));

    	// Stage 22 : Group by allele
    	BasicDBObject groupAllele = new BasicDBObject();
    	BasicDBObject groupAlleleId = new BasicDBObject();
    	groupAlleleId.put(FST_S22_VARIANTID, "$_id");
    	groupAlleleId.put(FST_S22_POPULATIONID, "$" + FST_S15_POPULATION);
    	groupAlleleId.put(FST_RES_ALLELEID, "$" + FST_S19_GENOTYPE);
    	groupAllele.put("_id", groupAlleleId);
    	groupAllele.put(FST_RES_SAMPLESIZE, new BasicDBObject("$first", "$" + FST_RES_SAMPLESIZE));
    	groupAllele.put(FST_S22_ALLELECOUNT, new BasicDBObject("$sum", 1));
    	groupAllele.put(FST_S22_HETEROZYGOTECOUNT, new BasicDBObject("$sum", new BasicDBObject("$toInt", "$" + FST_S20_HETEROZYGOTE)));
    	pipeline.add(new BasicDBObject("$group", groupAllele));

    	// Stage 23 : Group by population
    	BasicDBObject groupPopulation = new BasicDBObject();
    	BasicDBObject groupPopulationId = new BasicDBObject();
    	groupPopulationId.put(FST_S22_VARIANTID, "$_id." + FST_S22_VARIANTID);
    	groupPopulationId.put(FST_S22_POPULATIONID, "$_id." + FST_S22_POPULATIONID);
    	groupPopulation.put("_id", groupPopulationId);
    	groupPopulation.put(FST_RES_SAMPLESIZE, new BasicDBObject("$first", "$" + FST_RES_SAMPLESIZE));
    	groupPopulation.put(FST_RES_ALLELEMAX, new BasicDBObject("$max", "$_id." + FST_RES_ALLELEID));

    	BasicDBObject groupPopulationAllele = new BasicDBObject();
    	groupPopulationAllele.put(FST_RES_ALLELEID, "$_id." + FST_RES_ALLELEID);
    	BasicDBObject alleleFrequencyOperation = new BasicDBObject("$divide", Arrays.asList("$" + FST_S22_ALLELECOUNT, new BasicDBObject("$multiply", Arrays.asList("$" + FST_RES_SAMPLESIZE, 2))));
    	BasicDBObject hetFrequencyOperation = new BasicDBObject("$divide", Arrays.asList("$" + FST_S22_HETEROZYGOTECOUNT, "$" + FST_RES_SAMPLESIZE));
    	groupPopulationAllele.put(FST_RES_ALLELEFREQUENCY, alleleFrequencyOperation);
    	groupPopulationAllele.put(FST_RES_HETEROZYGOTEFREQUENCY, hetFrequencyOperation);

    	groupPopulation.put(FST_RES_ALLELES, new BasicDBObject("$push", groupPopulationAllele));
    	pipeline.add(new BasicDBObject("$group", groupPopulation));

    	// Stage 24 : Group by variant
    	BasicDBObject groupVariant = new BasicDBObject();
    	groupVariant.put("_id", "$_id." + FST_S22_VARIANTID);
    	groupVariant.put(FST_RES_ALLELEMAX, new BasicDBObject("$max", "$" + FST_RES_ALLELEMAX));
    	BasicDBObject groupVariantPopulation = new BasicDBObject();
    	groupVariantPopulation.put(FST_RES_SAMPLESIZE, "$" + FST_RES_SAMPLESIZE);
    	groupVariantPopulation.put(FST_RES_ALLELES, "$" + FST_RES_ALLELES);
    	groupVariant.put(FST_RES_POPULATIONS, new BasicDBObject("$push", groupVariantPopulation));
    	pipeline.add(new BasicDBObject("$group", groupVariant));

    	return pipeline;
    }

    private static final String TJD_S14_GENOTYPES = "gl";
    private static final String TJD_S15_SAMPLESIZE = "sz";
    private static final String TJD_S18_GENOTYPE = "gt";
    private static final String TJD_S20_VARIANTID = "vi";
    private static final String TJD_S20_ALLELEID = "ai";
    private static final String TJD_S20_ALLELECOUNT = "ac";
    private static final String TJD_S21_NUMALLELES = "na";
    private static final String TJD_S23_ALLELEFREQUENCY = "k";
    private static final String TJD_S24_FREQUENCYSUM = "kc";

    private static final String TJD_RES_SEGREGATINGSITES = "sg";
    private static final String TJD_RES_TAJIMAD = "tjd";

    private List<BasicDBObject> buildTajimaDQuery(GigwaDensityRequest gdr, boolean useTempColl) {
    	String info[] = GigwaSearchVariantsRequest.getInfoFromId(gdr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

    	List<String> selectedIndividuals = new ArrayList<String>();
        selectedIndividuals.addAll(gdr.getPlotIndividuals().size() == 0 ? MgdbDao.getProjectIndividuals(sModule, projId) : gdr.getPlotIndividuals().stream().map(csi -> csi.substring(1 + csi.lastIndexOf(IGigwaService.ID_SEPARATOR))).collect(Collectors.toSet()));

        TreeMap<String, List<GenotypingSample>> individualToSampleListMap = new TreeMap<String, List<GenotypingSample>>();
        individualToSampleListMap.putAll(MgdbDao.getSamplesByIndividualForProject(sModule, projId, selectedIndividuals));

        final int sampleSize = 2*selectedIndividuals.size();
        int intervalSize = Math.max(1, (int) ((gdr.getDisplayedRangeMax() - gdr.getDisplayedRangeMin()) / gdr.getDisplayedRangeIntervalCount()));
        List<Long> intervalBoundaries = new ArrayList<Long>();
        for (int i = 0; i < gdr.getDisplayedRangeIntervalCount(); i++) {
			intervalBoundaries.add(gdr.getDisplayedRangeMin() + (i*intervalSize));
		}
        intervalBoundaries.add(gdr.getDisplayedRangeMax() + 1);

        double a1 = 0, a2 = 0;
        for (int i = 1; i < sampleSize; i++) {
        	a1 += 1.0 / i;
        	a2 += 1.0 / (i*i);
        }

        double b1 = (double)(sampleSize + 1) / (double)(3*(sampleSize - 1));
        double b2 = 2.0*(sampleSize*sampleSize + sampleSize + 3) / (9.0*sampleSize*(sampleSize - 1));
        double c1 = b1 - 1/a1;
        double c2 = b2 - (double)(sampleSize + 2)/(a1*sampleSize) + a2/(a1*a1);
        double e1 = c1 / a1;
        double e2 = c2 / (a1*a1 + a2);

    	List<BasicDBObject> pipeline = buildGenotypeDataQuery(gdr, useTempColl, individualToSampleListMap, true);

    	// Stage 14 : Get the genotypes needed
    	BasicDBList genotypePaths = getFullPathToGenotypes(sModule, projId, selectedIndividuals, individualToSampleListMap);
    	BasicDBObject genotypeProjection = new BasicDBObject();
    	genotypeProjection.put(VariantData.FIELDNAME_REFERENCE_POSITION, 1);
    	genotypeProjection.put(TJD_S14_GENOTYPES, genotypePaths);
    	pipeline.add(new BasicDBObject("$project", genotypeProjection));

    	// Stage 15 : Count non-null genotypes
    	BasicDBObject sampleSizeMapping = new BasicDBObject();
    	sampleSizeMapping.put("input", "$" + TJD_S14_GENOTYPES);
    	sampleSizeMapping.put("in", new BasicDBObject("$cmp", Arrays.asList("$$this", null)));
    	BasicDBObject addSampleSize = new BasicDBObject("$sum", new BasicDBObject("$map", sampleSizeMapping));
    	pipeline.add(new BasicDBObject("$addFields", new BasicDBObject(TJD_S15_SAMPLESIZE, addSampleSize)));

    	// Stage 16 : Unwind individuals
    	pipeline.add(new BasicDBObject("$unwind", "$" + TJD_S14_GENOTYPES));

    	// Stage 17 : Eliminate missing genotypes
    	BasicDBObject matchMissing = new BasicDBObject(TJD_S14_GENOTYPES, new BasicDBObject("$ne", null));
    	pipeline.add(new BasicDBObject("$match", matchMissing));

    	// Stage 18 : Split the genotype string
    	BasicDBObject splitMapping = new BasicDBObject();
    	splitMapping.put("input", new BasicDBObject("$split", Arrays.asList("$" + TJD_S14_GENOTYPES, "/")));
    	splitMapping.put("in", new BasicDBObject("$toInt", "$$this"));
    	BasicDBObject splitProjection = new BasicDBObject();
    	splitProjection.put(VariantData.FIELDNAME_REFERENCE_POSITION, 1);
    	splitProjection.put(TJD_S18_GENOTYPE, new BasicDBObject("$map", splitMapping));
    	splitProjection.put(TJD_S15_SAMPLESIZE, 1);
    	pipeline.add(new BasicDBObject("$project", splitProjection));

    	// Stage 19 : Unwind alleles
    	pipeline.add(new BasicDBObject("$unwind", "$" + TJD_S18_GENOTYPE));

    	// Stage 20 : Count the alleles
    	BasicDBObject alleleGroupId = new BasicDBObject();
    	alleleGroupId.put(TJD_S20_VARIANTID, "$_id");
    	alleleGroupId.put(TJD_S20_ALLELEID, "$" + TJD_S18_GENOTYPE);
    	BasicDBObject alleleGroup = new BasicDBObject();
    	alleleGroup.put("_id", alleleGroupId);
    	alleleGroup.put(VariantData.FIELDNAME_REFERENCE_POSITION, new BasicDBObject("$first", "$" + VariantData.FIELDNAME_REFERENCE_POSITION));
    	alleleGroup.put(TJD_S20_ALLELECOUNT, new BasicDBObject("$sum", 1));
    	alleleGroup.put(TJD_S15_SAMPLESIZE, new BasicDBObject("$first", "$" + TJD_S15_SAMPLESIZE));
    	pipeline.add(new BasicDBObject("$group", alleleGroup));

    	// Stage 21 : Group by variant, keeping only one of the two alleles
    	BasicDBObject variantGroup = new BasicDBObject();
    	variantGroup.put("_id", "$_id." + TJD_S20_VARIANTID);
    	variantGroup.put(TJD_S20_ALLELECOUNT, new BasicDBObject("$first", "$" + TJD_S20_ALLELECOUNT));
    	variantGroup.put(TJD_S21_NUMALLELES, new BasicDBObject("$sum", 1));
    	variantGroup.put(TJD_S15_SAMPLESIZE, new BasicDBObject("$first", "$" + TJD_S15_SAMPLESIZE));
    	variantGroup.put(VariantData.FIELDNAME_REFERENCE_POSITION, new BasicDBObject("$first", "$" + VariantData.FIELDNAME_REFERENCE_POSITION));
    	pipeline.add(new BasicDBObject("$group", variantGroup));

    	// Stage 22 : Keep only biallelic variants
    	pipeline.add(new BasicDBObject("$match", new BasicDBObject(TJD_S21_NUMALLELES, 2)));

    	// Stage 23 : Compute the average pairwise polymorphism for one variant
    	BasicDBObject alleleFrequency = new BasicDBObject();
    	alleleFrequency.put("vars", new BasicDBObject("freq",
    			new BasicDBObject("$divide", Arrays.asList("$" + TJD_S20_ALLELECOUNT,
    				new BasicDBObject("$multiply", Arrays.asList("$" + TJD_S15_SAMPLESIZE, 2))))));
    	alleleFrequency.put("in", new BasicDBObject("$multiply", Arrays.asList("$$freq", new BasicDBObject("$subtract", Arrays.asList(1, "$$freq")))));  // p(1-p)
    	BasicDBObject frequencyProjection = new BasicDBObject();
    	frequencyProjection.put(VariantData.FIELDNAME_REFERENCE_POSITION, 1);
    	frequencyProjection.put(TJD_S23_ALLELEFREQUENCY, new BasicDBObject("$let", alleleFrequency));
    	pipeline.add(new BasicDBObject("$project", frequencyProjection));

    	// Stage 24 : Group by graph interval
    	BasicDBObject outputGroup = new BasicDBObject();
    	outputGroup.put("_id", 1);
    	outputGroup.put(TJD_S24_FREQUENCYSUM, new BasicDBObject("$sum", "$" + TJD_S23_ALLELEFREQUENCY));
    	outputGroup.put(TJD_RES_SEGREGATINGSITES, new BasicDBObject("$sum", 1));
    	pipeline.add(new BasicDBObject("$group", outputGroup));

    	// Stage 25 : Compute the Tajima's D value
    	BasicDBObject finalProject = new BasicDBObject();
    	finalProject.put(TJD_RES_SEGREGATINGSITES, 1);
    	finalProject.put(TJD_RES_TAJIMAD, new BasicDBObject("$divide", Arrays.asList(
    		new BasicDBObject("$subtract", Arrays.asList(
    			new BasicDBObject("$divide", Arrays.asList(
    				new BasicDBObject("$multiply", Arrays.asList("$" + TJD_S24_FREQUENCYSUM, 2*sampleSize)),
    				sampleSize - 1
    			)),
    			new BasicDBObject("$divide", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, a1))
    		)),
    		new BasicDBObject("$sqrt", new BasicDBObject("$abs", new BasicDBObject("$add", Arrays.asList(
    			new BasicDBObject("$multiply", Arrays.asList(e1, "$" + TJD_RES_SEGREGATINGSITES)),
    			new BasicDBObject("$multiply", Arrays.asList(
    				new BasicDBObject("$multiply", Arrays.asList(e2, "$" + TJD_RES_SEGREGATINGSITES)),
    				new BasicDBObject("$subtract", Arrays.asList("$" + TJD_RES_SEGREGATINGSITES, 1))
    			))
    		))))
    	)));
    	pipeline.add(new BasicDBObject("$project", finalProject));

    	return pipeline;
    }


    private BasicDBList getFullPathToGenotypes(String sModule, int projId, Collection<String> selectedIndividuals, Map<String, List<GenotypingSample>> individualToSampleListMap){
    	BasicDBList result = new BasicDBList();
    	Iterator<String> indIt = selectedIndividuals.iterator();
        while (indIt.hasNext()) {
            String individual = indIt.next();
            List<GenotypingSample> individualSamples = individualToSampleListMap.get(individual);

            int finalSample = Integer.MAX_VALUE;
            for (int k=0; k<individualSamples.size(); k++) {    // this loop is executed only once for single-run projects
                GenotypingSample individualSample = individualSamples.get(k);
                if (individualSample.getId() < finalSample)
                	finalSample = individualSample.getId();
            }

            String pathToGT = finalSample + "." + SampleGenotype.FIELDNAME_GENOTYPECODE;
            String fullPathToGT = "$" + VariantRunData.FIELDNAME_SAMPLEGENOTYPES/* + (int) ((individualSample.getId() - 1) / 100)*/ + "." + pathToGT;
            result.add(fullPathToGT);
        }
        return result;
    }

    
    public Map<Long, Integer> selectionVcfFieldPlotData(GigwaVcfFieldPlotRequest gvfpr) throws Exception {
        long before = System.currentTimeMillis();

        String info[] = GigwaSearchVariantsRequest.getInfoFromId(gvfpr.getVariantSetId(), 2);
        String sModule = info[0];
        int projId = Integer.parseInt(info[1]);

        ProgressIndicator progress = new ProgressIndicator(tokenManager.readToken(gvfpr.getRequest()), new String[] {"Calculating plot data for " + gvfpr.getVcfField() +  " field regarding " + (gvfpr.getDisplayedVariantType() != null ? gvfpr.getDisplayedVariantType() + " " : "") + "variants on sequence " + gvfpr.getDisplayedSequence()});
        ProgressIndicator.registerProgressIndicator(progress);

        final MongoTemplate mongoTemplate = MongoTemplateManager.get(sModule);
        Collection<BasicDBList> variantQueryDBListColl = ga4ghService.buildVariantDataQuery(gvfpr, ga4ghService.getSequenceIDsBeingFilteredOn(gvfpr.getRequest().getSession(), sModule), true);
        final BasicDBList variantQueryDBList = variantQueryDBListColl.size() == 1 ? variantQueryDBListColl.iterator().next() : new BasicDBList();

        MongoCollection<Document> tmpVarColl = ga4ghService.getTemporaryVariantCollection(sModule, tokenManager.readToken(gvfpr.getRequest()), false);
        long nTempVarCount = mongoTemplate.count(new Query(), tmpVarColl.getNamespace().getCollectionName());
        if (GenotypingDataQueryBuilder.getGroupsForWhichToFilterOnGenotypingOrAnnotationData(gvfpr, false).size() > 0 && nTempVarCount == 0)
        {
            progress.setError(Constants.MESSAGE_TEMP_RECORDS_NOT_FOUND);
            return null;
        }

        final String mainVarCollName = mongoTemplate.getCollectionName(VariantData.class), usedVarCollName = nTempVarCount == 0 ? mainVarCollName : tmpVarColl.getNamespace().getCollectionName();
        final ConcurrentHashMap<Long, Integer> result = new ConcurrentHashMap<Long, Integer>();

        if (gvfpr.getDisplayedRangeMin() == null || gvfpr.getDisplayedRangeMax() == null)
            if (!findDefaultRangeMinMax(gvfpr, usedVarCollName, progress))
                return result;

		final int intervalSize = Math.max(1, (int) ((gvfpr.getDisplayedRangeMax() - gvfpr.getDisplayedRangeMin()) / gvfpr.getDisplayedRangeIntervalCount()));
		final ArrayList<Thread> threadsToWaitFor = new ArrayList<Thread>();
		final long rangeMin = gvfpr.getDisplayedRangeMin();
		final ProgressIndicator finalProgress = progress;

		if (gvfpr.getPlotIndividuals() == null || gvfpr.getPlotIndividuals().size() == 0)
			gvfpr.setPlotIndividuals(MgdbDao.getProjectIndividuals(sModule, projId));

        List<Integer>[] sampleIDsGroupedBySortedIndividuals = new List[gvfpr.getPlotIndividuals().size()];
        TreeMap<String, ArrayList<GenotypingSample>> samplesByIndividual = MgdbDao.getSamplesByIndividualForProject(sModule, projId, gvfpr.getPlotIndividuals());
        int k = 0;
        for (String ind : gvfpr.getPlotIndividuals()) {
        	sampleIDsGroupedBySortedIndividuals[k] = samplesByIndividual.get(ind).stream().map(sp -> sp.getId()).collect(Collectors.toList());
            k++;
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
                        List<String> variantsInInterval = mongoTemplate.getCollection(usedVarCollName).distinct("_id", query.getQueryObject(), String.class).into(new ArrayList<>());

                        final ArrayList<BasicDBObject> pipeline = new ArrayList<BasicDBObject>();

            			BasicDBObject group = new BasicDBObject();
            			ArrayList<Object> individualValuePaths = new ArrayList<>();
            			for (int j=0; j<gvfpr.getPlotIndividuals().size(); j++)
            			{
            				List<Integer> individualSamples = sampleIDsGroupedBySortedIndividuals[j];
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

	        			Iterator<Document> it = mongoTemplate.getCollection(mongoTemplate.getCollectionName(VariantRunData.class)).aggregate(pipeline).iterator();
	        			result.put(rangeMin + (chunkIndex*intervalSize), it.hasNext() ? Double.valueOf(it.next().get(gvfpr.getVcfField()).toString()).intValue() : 0);
	        			finalProgress.setCurrentStepProgress((short) result.size() * 100 / gvfpr.getDisplayedRangeIntervalCount());
            		}
            	}
            };

            if (chunkIndex%GigwaGa4ghServiceImpl.INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS  == (GigwaGa4ghServiceImpl.INITIAL_NUMBER_OF_SIMULTANEOUS_QUERY_THREADS-1))
                t.run();    // run synchronously
            else
            {
                threadsToWaitFor.add(t);
                t.start();    // run asynchronously for better speed
            }
        }

        if (progress.isAborted())
            return null;

        for (Thread ttwf : threadsToWaitFor)    // wait for all threads before moving to next phase
            ttwf.join();

        progress.setCurrentStepProgress(100);
        LOG.debug("selectionVcfFieldPlotData treated " + gvfpr.getDisplayedRangeIntervalCount() + " intervals on sequence " + gvfpr.getDisplayedSequence() + " between " + gvfpr.getDisplayedRangeMin() + " and " + gvfpr.getDisplayedRangeMax() + " in " + (System.currentTimeMillis() - before)/1000f + "s");
        progress.markAsComplete();

        return new TreeMap<Long, Integer>(result);
    }
}
