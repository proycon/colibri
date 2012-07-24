
class LanguageModel {
    private:
        int order;
    public:
        std::unordered_map<EncNGram> ngrams;
        std::unordered_map<EncNGram> backoff;
        std::map<int,unsigned int> total;
        
        LanguageModel(const std::string & filename,  ClassEncoder & encoder);
                
        
        double score(EncData & data, bool fullsentence = false);
        double score(EncNGram & ngram);
};
