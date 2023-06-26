#define PROFILE
#include "openfhe.h"

using namespace lbcrypto;

using namespace std;


void innerProduct(Ciphertext<DCRTPoly>& output, const Ciphertext<DCRTPoly>& input1, const Plaintext& input2, CryptoContext<DCRTPoly>& cc, int i, bool accumulate=true){
	
    auto cMul = cc->EvalMult(input1, input2);
    auto ctemp2= cMul;
	
    for(int i=1; i<=log(8)+1; i++)
    {
      	auto ctemp1 = cc->EvalRotate(cMul, i);
	ctemp2 = cc->EvalAdd(ctemp2, ctemp1);
    }
    auto ctemp1=cc->EvalRotate(ctemp2, 4);
    ctemp2=cc->EvalAdd(ctemp2, ctemp1);

    // apply the mask
    std::vector<double> x1 = {0, 0, 0, 0, 0, 0,0, 0};
    x1[i]=1;
    
    // put the result in output
    Plaintext mask = cc->MakeCKKSPackedPlaintext(x1);
    if (accumulate){
    	output += cc->EvalMult(ctemp2, mask);
    }else{
    	output = cc->EvalMult(ctemp2, mask);
    }
}


void EvalLogisticExample(Ciphertext<DCRTPoly>& output, const Ciphertext<DCRTPoly>& input1,  CryptoContext<DCRTPoly>& cc, double a, double b) {

    uint32_t polyDegree = 3;
    double lowerBound = -a;
    double upperBound = b;
    
    //output  = cc->EvalLogistic(input1, lowerBound, upperBound, polyDegree);
    output = cc->EvalChebyshevFunction([](double x) -> double { if (x < 0) return 0; else return x; }, input1, lowerBound,
                                            upperBound, polyDegree);
}


float randomFloat()
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    return dis(gen);
}


//  --- BEGIN   functions that print keys  -----------------------------------------//
void print_moduli_chain(const DCRTPoly& poly){
    int num_primes = poly.GetNumOfElements();
    double total_bit_len = 0.0;
    for (int i = 0; i < num_primes; i++) {
        auto qi = poly.GetParams()->GetParams()[i]->GetModulus();
        std::cout << "q_" << i << ": " 
                    << qi
                    << ",  log q_" << i <<": " << log(qi.ConvertToDouble()) / log(2)
                    << std::endl;
        total_bit_len += log(qi.ConvertToDouble()) / log(2);
    }
    std::cout << "Total bit length: " << total_bit_len << std::endl;
}

void print_in_coeff_domain(const DCRTPoly& poly_rns){
    std::cout << poly_rns.CRTInterpolate() << std::endl;
}

/* XXX: assuming that sk has coefficients in {-1, 0, 1} */
void print_ternary_sk(const KeyPair<DCRTPoly>& keys){
    const DCRTPoly& ckks_sk = keys.secretKey->GetPrivateElement();

    Poly poly = ckks_sk.CRTInterpolate();
    const BigInteger Q = poly.GetModulus();
    const BigInteger Qover2 = Q / 2;
    int N = poly.GetRingDimension();
    std::cout << "[ ";
    for(int i = 0; i < N; i++){
        const BigInteger& ai = poly.at(i);
        if(ai > Qover2)
            std::cout << "-1 ";
        else
            std::cout << ai << " ";
    }
    std::cout << "]" << std::endl;
}
void print_rlwe_sample_in_coeff_domain(const std::vector<DCRTPoly>& rlwe_sample){
    std::cout << "  -a: " << rlwe_sample[1].CRTInterpolate() << std::endl;
    std::cout << "   b: " << rlwe_sample[0].CRTInterpolate() << std::endl;
}

void print_rlwe_sample_in_rns_domain(const std::vector<DCRTPoly>& rlwe_sample){
    std::cout << "  -a: " << rlwe_sample[1] << std::endl;
    std::cout << "   b: " << rlwe_sample[0] << std::endl;
}

void print_relin_key(const EvalKey<DCRTPoly>& relinKey){
    const std::vector<DCRTPoly>& relinKey_row_a = relinKey->GetAVector();
    std::cout << "Modulus chain of relinearization key:" << std::endl;
    print_moduli_chain(relinKey_row_a[0]);

    #if VERBOSE
    const std::vector<DCRTPoly>& relinKey_row_b = relinKey->GetBVector();
    int nrows_ks = relinKey_row_a.size();
    for(int i = 0; i < nrows_ks; i++){
        std::cout << "row " << i << " of relin key:" << std::endl;
        std::cout << "    a_" << i << ": " << relinKey_row_a[i] << std::endl;
        std::cout << "    b_" << i << ": " << relinKey_row_b[i] << std::endl;
    }
    #endif
}

void print_automorphism_key(const EvalKey<DCRTPoly>& autKey){
    const std::vector<DCRTPoly>& row_b = autKey->GetBVector();
    std::cout << "Modulus chain of automorphism key:" << std::endl;
    print_moduli_chain(row_b[0]);

    #if VERBOSE
    const std::vector<DCRTPoly>& row_a = autKey->GetAVector();
    int nrows_ks = row_a.size();
    for(int i = 0; i < nrows_ks; i++){
        std::cout << "row " << i << " of automorphism key:" << std::endl;
        std::cout << "    a_" << i << ": " << row_a[i] << std::endl;
        std::cout << "    b_" << i << ": " << row_b[i] << std::endl;
    }
    #endif
}

void print_all_automorphism_keys(const std::map<usint, EvalKey<DCRTPoly>>& aut_key_map){
    std::cout << "# of automorphism keys: " << aut_key_map.size() << std::endl;
    for(auto const &ent : aut_key_map) {
        int k = ent.first; // int k corresponding to automorphism X --> X^k
        const EvalKey<DCRTPoly>& evk_i = ent.second; // key to be used after automorphism X --> X^k
        std::cout << "automorphism key for X --> X^"<< k << std::endl;
        print_automorphism_key(evk_i);
        std::cout << std::endl;
    }
}
//  --- END of functions that print keys  -----------------------------------------//


int main() {

    int n;
    cout << "enter the number of hidden layer" << endl;
    cin >> n ;
    cout<<"This is a "<<n+1<<" layer Neural Metwork"<<endl;
    
    
    uint32_t multDepth = (n+1)*5+2;
    uint32_t scaleModSize =24;
    uint32_t batchSize = 8;
    uint32_t firstModSize   = 25;
    
    
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);
    parameters.SetBatchSize(batchSize);
    parameters.SetFirstModSize(firstModSize);
    parameters.SetSecurityLevel(HEStd_128_classic);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    parameters.SetMultiplicativeDepth(multDepth);


    // Enable the features that you wish to use
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);


    // B. Step 2: Key Generation
    cout << "keyGeneration starts" << endl;
    auto keyPair = cc->KeyGen();
    cout << "keyGeneration ends" << endl;
    
    
    cout << "Relin keyGeneration starts" << endl;
    cc->EvalMultKeyGen(keyPair.secretKey);
    cout << "Relin keyGeneration ends" << endl;
    
    
    cout << "Rotation keyGeneration starts" << endl;
    cc->EvalRotateKeyGen(keyPair.secretKey, {1, 2,3, 4, 5, 6, 8, 7});
    cout << "Rotation keyGeneration ends" << endl;
    
    
    cout << "Encryption of the plaintext starts" << endl;
    std::vector<double> x1 = {0.1, 0.5, 0.75, 0.1, 0.12, 0.3, 0.4, 0.25};
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);
    auto c1 = cc->Encrypt(keyPair.publicKey, ptxt1);
    cout << "Encryption of the plaintext ends" << endl;
    
    
    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << std::endl << std::endl;



    //print relin key
    const std::string& id_sk = keyPair.secretKey->GetKeyTag();
    const std::vector<EvalKey<DCRTPoly>>& relinKeys = cc->GetEvalMultKeyVector(id_sk);
    std::cout << "# of relin keys: " << relinKeys.size() << std::endl;
    if (relinKeys.size() > 1){
        std::cout << "ERROR: more than one relinearization key created." << std::endl;
        exit(1);
    }
    const EvalKey<DCRTPoly>& relinKey = relinKeys[0];
    print_relin_key(relinKey);
    
    
    Ciphertext<DCRTPoly> output, temp, output_f;
    std::vector<double> w={0, 0, 0, 0, 0, 0, 0, 0};
    for(int k=0; k<8; k++)
    {
    	
    	w[k]=randomFloat();
    }
    Plaintext pw = cc->MakeCKKSPackedPlaintext(w);
    
    
    std::cout << "Ptxt precision before computation starts: " << ptxt1->GetLogPrecision() << std::endl;
    
    
    std::cout << "Number of levels remaining at the start of input layer: " << multDepth - c1->GetLevel() << std::endl;
    cout<<"Comptation for the input layer starts"<<endl;
    innerProduct(temp, c1, pw, cc, 0, false);
    innerProduct(temp, c1, pw, cc, 1);
    innerProduct(temp, c1, pw, cc, 2);
    innerProduct(temp, c1, pw, cc, 3);
    innerProduct(temp, c1, pw, cc, 4);
    innerProduct(temp, c1, pw, cc, 5);
    innerProduct(temp, c1, pw, cc, 6);
    innerProduct(temp, c1, pw, cc, 7);
    std::cout << "Number of levels remaining after th innerproduct: " << multDepth - temp->GetLevel() << std::endl;
    EvalLogisticExample(output, temp, cc, sqrt(8), sqrt(8));
    cout<<"Comptation for the input layer ends"<<endl;
    std::cout << "Number of levels remaining at the end of input layer: " << multDepth - output->GetLevel() << std::endl;
    
    
    for(int i =0; i < n; i++)
    {
    	std::cout << "Number of levels remaining at the begining of the "<<i+2<<"th layer " << multDepth - output->GetLevel() << std::endl;
    	cout<<"Comptation for the "<<i+2<<"th layer starts"<<endl;
    	
    	
    	for(int k=0; k<8; k++)
    	{
    		w[k]=randomFloat();
    	}
    	pw = cc->MakeCKKSPackedPlaintext(w);
    	innerProduct(temp, output, pw, cc, 0);

    	
    	for(int k=0; k<8; k++)
    	{
    		w[k]=randomFloat();
    	}
    	pw = cc->MakeCKKSPackedPlaintext(w);
    	innerProduct(temp, output, pw, cc, 1);
    	
    	
    	for(int k=0; k<8; k++)
    	{
    		w[k]=randomFloat();
    	}
    	pw = cc->MakeCKKSPackedPlaintext(w);
    	innerProduct(temp, output, pw, cc, 2);
    	
    	
    	for(int k=0; k<8; k++)
    	{
    		w[k]=randomFloat();
    	}
    	pw = cc->MakeCKKSPackedPlaintext(w);
    	innerProduct(temp, output, pw, cc, 3);
    	
    	
    	for(int k=0; k<8; k++)
    	{
    		w[k]=randomFloat();
    	}
    	pw = cc->MakeCKKSPackedPlaintext(w);
    	innerProduct(temp, output, pw, cc, 4);
    	
    	
    	for(int k=0; k<8; k++)
    	{
    		w[k]=randomFloat();
    	}
    	pw = cc->MakeCKKSPackedPlaintext(w);
    	innerProduct(temp, output, pw, cc, 5);
    	
    	
    	for(int k=0; k<8; k++)
    	{
    		w[k]=randomFloat();
    	}
    	pw = cc->MakeCKKSPackedPlaintext(w);
    	innerProduct(temp, output, pw, cc, 6);
    	
    	
    	for(int k=0; k<8; k++)
    	{
    		w[k]=randomFloat();
    	}
    	pw = cc->MakeCKKSPackedPlaintext(w);
    	innerProduct(temp, output, pw, cc, 7);

	
    	EvalLogisticExample(output, temp, cc, pow(sqrt(8),i+2), pow(sqrt(8), i+2));
 
    	
    	cout<<"Comptation for the "<<i+2<<"th layer ends"<<endl;
    	std::cout << "Number of levels remaining at the end of the "<<i+2<<"th layer " << multDepth - output->GetLevel() << std::endl;
    }


    std::vector<double> w16 = {0.25, 0.5, 0.75, 0.1, 0.2, 0.3, 0.4, 0.5};
    Plaintext w_f = cc->MakeCKKSPackedPlaintext(w16);
    
    
    std::cout << "Number of levels remaining at the start of output layer: " << multDepth - output->GetLevel() << std::endl;
    cout<<"Comptation for the output layer starts"<<endl;
    innerProduct(output_f, output, w_f, cc, 0, false);
    cout<<"Comptation for the output layer ends"<<endl;
    std::cout << "Number of levels remaining at the end of output layer: " << multDepth - output_f->GetLevel() << std::endl;
    
    
    cout<< "Decryption starts here " << endl;
    Plaintext plaintextDec;
    cc->Decrypt(keyPair.secretKey, output_f, &plaintextDec);
    plaintextDec->SetLength(8);
    
    
    std::vector<std::complex<double>> finalResult = plaintextDec->GetCKKSPackedValue();
    std::cout << "Final output\n\t" << finalResult[0] << std::endl << std::endl;
    std::cout << "Ptxt precision after decryption: " << plaintextDec->GetLogPrecision() << std::endl;
    
        
    return 0;
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
