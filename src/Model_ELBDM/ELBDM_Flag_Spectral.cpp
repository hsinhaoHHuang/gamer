#include "GAMER.h"

#if ( MODEL == ELBDM )

#include <complex.h>

#define flag_spectral_float double
typedef std::complex<real>                complex_type;
typedef std::complex<flag_spectral_float> flag_spectral_complex_type;


// GramFE parameters for extension of size 14
#define FLAG_SPECTRAL_ORDER   14
#define FLAG_SPECTRAL_NDELTA  14
#define FLAG_SPECTRAL_ND      32

const static flag_spectral_float  Flag_Spectral_Polynomials[FLAG_SPECTRAL_ORDER][FLAG_SPECTRAL_NDELTA] = {
{0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000, 0.2672612419124243965384835064469370990991592407226562500000000000},
{-0.4309458036856673168735198942158604040741920471191406250000000000, -0.3646464492724877382023862537607783451676368713378906250000000000, -0.2983470948593081595312526133056962862610816955566406250000000000, -0.2320477404461285531045433572217007167637348175048828125000000000, -0.1657483860329489744334097167666186578571796417236328125000000000, -0.0994490316197693818844882684970798436552286148071289062500000000, -0.0331496772065897962744607241347694071009755134582519531250000000, 0.0331496772065897962744607241347694071009755134582519531250000000, 0.0994490316197693818844882684970798436552286148071289062500000000, 0.1657483860329489744334097167666186578571796417236328125000000000, 0.2320477404461285531045433572217007167637348175048828125000000000, 0.2983470948593081595312526133056962862610816955566406250000000000, 0.3646464492724877382023862537607783451676368713378906250000000000, 0.4309458036856673168735198942158604040741920471191406250000000000},
{0.4818120558297157574045854744326788932085037231445312500000000000, 0.2594372608313854078332383323868270963430404663085937500000000000, 0.0741249316661101165237823806819505989551544189453125000000000000, -0.0741249316661101165237823806819505989551544189453125000000000000, -0.1853123291652752913094559517048764973878860473632812500000000000, -0.2594372608313854078332383323868270963430404663085937500000000000, -0.2964997266644404660951295227278023958206176757812500000000000000, -0.2964997266644404660951295227278023958206176757812500000000000000, -0.2594372608313854078332383323868270963430404663085937500000000000, -0.1853123291652752913094559517048764973878860473632812500000000000, -0.0741249316661101165237823806819505989551544189453125000000000000, 0.0741249316661101165237823806819505989551544189453125000000000000, 0.2594372608313854078332383323868270963430404663085937500000000000, 0.4818120558297157574045854744326788932085037231445312500000000000},
{-0.4585783658733355583336788185988552868366241455078125000000000000, -0.0352752589133335028859228543751669349148869514465332031250000000, 0.2116515534800010311933249340654583647847175598144531250000000000, 0.3142704885006075699038774473592638969421386718750000000000000000, 0.3046499633424257225122744330292334780097007751464843750000000000, 0.2148583951993949803238592721754685044288635253906250000000000000, 0.0769642012654549179107021927848109044134616851806640625000000000, -0.0769642012654549179107021927848109044134616851806640625000000000, -0.2148583951993949803238592721754685044288635253906250000000000000, -0.3046499633424257225122744330292334780097007751464843750000000000, -0.3142704885006075699038774473592638969421386718750000000000000000, -0.2116515534800010311933249340654583647847175598144531250000000000, 0.0352752589133335028859228543751669349148869514465332031250000000, 0.4585783658733355583336788185988552868366241455078125000000000000},
{0.3875694570442999031811837085115257650613784790039062500000000000, -0.2086912461007768837539799733349354937672615051269531250000000000, -0.3577564218870460943655587016110075637698173522949218750000000000, -0.2493453849515775699874353676932514645159244537353515625000000000, -0.0352335870040272683412219123511022189632058143615722656250000000, 0.1707473831733629099360882719338405877351760864257812500000000000, 0.2927097997257649963920300706377020105719566345214843750000000000, 0.2927097997257649963920300706377020105719566345214843750000000000, 0.1707473831733629099360882719338405877351760864257812500000000000, -0.0352335870040272683412219123511022189632058143615722656250000000, -0.2493453849515775699874353676932514645159244537353515625000000000, -0.3577564218870460943655587016110075637698173522949218750000000000, -0.2086912461007768837539799733349354937672615051269531250000000000, 0.3875694570442999031811837085115257650613784790039062500000000000},
{-0.2948961391092899120280890201684087514877319335937500000000000000, 0.3856334126813791285393051566643407568335533142089843750000000000, 0.2722118207162676495336484094877960160374641418457031250000000000, -0.0577419013640567690970328840194270014762878417968750000000000000, -0.2866472960572818418079066304926527664065361022949218750000000000, -0.2990205606352939748937558306352002546191215515136718750000000000, -0.1237326457801216500476115811579802539199590682983398437500000000, 0.1237326457801216500476115811579802539199590682983398437500000000, 0.2990205606352939748937558306352002546191215515136718750000000000, 0.2866472960572818418079066304926527664065361022949218750000000000, 0.0577419013640567690970328840194270014762878417968750000000000000, -0.2722118207162676495336484094877960160374641418457031250000000000, -0.3856334126813791285393051566643407568335533142089843750000000000, 0.2948961391092899120280890201684087514877319335937500000000000000},
{0.2027563273040598468277551091887289658188819885253906250000000000, -0.4523025762936719873508195632894057780504226684570312500000000000, -0.0155966405618507587826915283812923007644712924957275390625000000, 0.3218579461400111196844875394162954762578010559082031250000000000, 0.2623071367220354832561213243025122210383415222167968750000000000, -0.0354469103678426294967707121941202785819768905639648437500000000, -0.2835752829427410359741656975529622286558151245117187500000000000, -0.2835752829427410359741656975529622286558151245117187500000000000, -0.0354469103678426294967707121941202785819768905639648437500000000, 0.2623071367220354832561213243025122210383415222167968750000000000, 0.3218579461400111196844875394162954762578010559082031250000000000, -0.0155966405618507587826915283812923007644712924957275390625000000, -0.4523025762936719873508195632894057780504226684570312500000000000, 0.2027563273040598468277551091887289658188819885253906250000000000},
{-0.1257441362172087295778766247167368419468402862548828125000000000, 0.4159229121030749709575502492953091859817504882812500000000000000, -0.2611608982972796755284150549414334818720817565917968750000000000, -0.3104033572354872316800822318327846005558967590332031250000000000, 0.0835363142701736199891016099172702524811029434204101562500000000, 0.3297486089612116644254058428487041965126991271972656250000000000, 0.1758659247793128987957800291042076423764228820800781250000000000, -0.1758659247793128987957800291042076423764228820800781250000000000, -0.3297486089612116644254058428487041965126991271972656250000000000, -0.0835363142701736199891016099172702524811029434204101562500000000, 0.3104033572354872316800822318327846005558967590332031250000000000, 0.2611608982972796755284150549414334818720817565917968750000000000, -0.4159229121030749709575502492953091859817504882812500000000000000, 0.1257441362172087295778766247167368419468402862548828125000000000},
{0.0699086407042275592704783093722653575241565704345703125000000000, -0.3172776770422635062018912321946118026971817016601562500000000000, 0.4248294319718443623479231519013410434126853942871093750000000000, 0.0376431142253532968755536103344638831913471221923828125000000000, -0.3495432035211377685968159312324132770299911499023437500000000000, -0.1344396936619760840603277074478683061897754669189453125000000000, 0.2688793873239521681206554148957366123795509338378906250000000000, 0.2688793873239521681206554148957366123795509338378906250000000000, -0.1344396936619760840603277074478683061897754669189453125000000000, -0.3495432035211377685968159312324132770299911499023437500000000000, 0.0376431142253532968755536103344638831913471221923828125000000000, 0.4248294319718443623479231519013410434126853942871093750000000000, -0.3172776770422635062018912321946118026971817016601562500000000000, 0.0699086407042275592704783093722653575241565704345703125000000000},
{-0.0344591278812576007339885109104216098785400390625000000000000000, 0.2041040651428334962158572807311429642140865325927734375000000000, -0.4320644495880761049022567021893337368965148925781250000000000000, 0.2836251294841972137028562883642734959721565246582031250000000000, 0.2359124908793789887617720069101778790354728698730468750000000000, -0.2783237251947729418155574876436730846762657165527343750000000000, -0.2385631930240910969498457916415645740926265716552734375000000000, 0.2385631930240910969498457916415645740926265716552734375000000000, 0.2783237251947729418155574876436730846762657165527343750000000000, -0.2359124908793789887617720069101778790354728698730468750000000000, -0.2836251294841972137028562883642734959721565246582031250000000000, 0.4320644495880761049022567021893337368965148925781250000000000000, -0.2041040651428334962158572807311429642140865325927734375000000000, 0.0344591278812576007339885109104216098785400390625000000000000000},
{0.0147897728358395794817647939112248423043638467788696289062500000, -0.1103544588520337932369130840015714056789875030517578125000000000, 0.3276503520555230086763742747280048206448554992675781250000000000, -0.4459685347422396306527048182033468037843704223632812500000000000, 0.1422093541907651914613097687833942472934722900390625000000000000, 0.3174112785537879233288549585267901420593261718750000000000000000, -0.2457377640416422426294928982315468601882457733154296875000000000, -0.2457377640416422426294928982315468601882457733154296875000000000, 0.3174112785537879233288549585267901420593261718750000000000000000, 0.1422093541907651914613097687833942472934722900390625000000000000, -0.4459685347422396306527048182033468037843704223632812500000000000, 0.3276503520555230086763742747280048206448554992675781250000000000, -0.1103544588520337932369130840015714056789875030517578125000000000, 0.0147897728358395794817647939112248423043638467788696289062500000},
{-0.0053617479838052716639706929413478064816445112228393554687500000, 0.0490806161594482537324779514165129512548446655273437500000000000, -0.1913731588065881450422267562316847033798694610595703125000000000, 0.3992440037171925415471207543305354192852973937988281250000000000, -0.4310020494674237645504888405412202700972557067871093750000000000, 0.0952741372506936690101042586320545524358749389648437500000000000, 0.3266541848595211350314571063790936022996902465820312500000000000, -0.3266541848595211350314571063790936022996902465820312500000000000, -0.0952741372506936690101042586320545524358749389648437500000000000, 0.4310020494674237645504888405412202700972557067871093750000000000, -0.3992440037171925415471207543305354192852973937988281250000000000, 0.1913731588065881450422267562316847033798694610595703125000000000, -0.0490806161594482537324779514165129512548446655273437500000000000, 0.0053617479838052716639706929413478064816445112228393554687500000},
{0.0015503894602372355823738381275234132772311568260192871093750000, -0.0170542840626095905387504814143539988435804843902587890625000000, 0.0837210308528107266523576868166856002062559127807617187500000000, -0.2387599768765342744814006437081843614578247070312500000000000000, 0.4263571015652398155104663146630628034472465515136718750000000000, -0.4604656696904589896490733735845424234867095947265625000000000000, 0.2046514087513151003427935847867047414183616638183593750000000000, 0.2046514087513151003427935847867047414183616638183593750000000000, -0.4604656696904589896490733735845424234867095947265625000000000000, 0.4263571015652398155104663146630628034472465515136718750000000000, -0.2387599768765342744814006437081843614578247070312500000000000000, 0.0837210308528107266523576868166856002062559127807617187500000000, -0.0170542840626095905387504814143539988435804843902587890625000000, 0.0015503894602372355823738381275234132772311568260192871093750000},
{-0.0003100778920474471056327459006496383153717033565044403076171875, 0.0040310125966168128611166743269222934031859040260314941406250000, -0.0241860755797008754319765699847266660071909427642822265625000000, 0.0886822771255698777403964072618691716343164443969726562500000000, -0.2217056928139246874120971142474445514380931854248046875000000000, 0.3990702470650644428928899287711828947067260742187500000000000000, -0.5320936627534192941979540592001285403966903686523437500000000000, 0.5320936627534192941979540592001285403966903686523437500000000000, -0.3990702470650644428928899287711828947067260742187500000000000000, 0.2217056928139246874120971142474445514380931854248046875000000000, -0.0886822771255698777403964072618691716343164443969726562500000000, 0.0241860755797008754319765699847266660071909427642822265625000000, -0.0040310125966168128611166743269222934031859040260314941406250000, 0.0003100778920474471056327459006496383153717033565044403076171875},
};

//-------------------------------------------------------------------------------------------------------
// Function    :  Least_Squares_Regression
// Description :  Compute the least squares linear regression for a set of data points (x, y).
//                The function calculates the slope (m) and intercept (b) of the linear regression model
//                using the least squares method.
//
// Parameter   :  x        : Array of x-coordinates of the data points
//                y        : Array of y-coordinates of the data points
//                n_start  : Starting index for linear regression
//                n_end    : End index for linear regression (set to size of x and y to include last point)
//                slope    : Pointer to the variable where the computed slope (m) will be stored
//                intercept: Pointer to the variable where the computed intercept (b) will be stored
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Least_Squares_Regression(flag_spectral_float x[], flag_spectral_float y[], int n_start, int n_end, flag_spectral_float *slope, flag_spectral_float *intercept) {
   flag_spectral_float sum_x = 0.0, sum_y = 0.0, sum_x_squared = 0.0, sum_xy = 0.0;

   const int n = n_end - n_start;

   for (int i = n_start; i < n_end; i++) {
      sum_x += x[i];
      sum_y += y[i];
      sum_x_squared += x[i] * x[i];
      sum_xy += x[i] * y[i];
   }

   *slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x * sum_x);
   *intercept = (sum_y - (*slope) * sum_x) / n;
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Prepare_for_Spectral_Criterion
// Description :  Evaluate decay of the coefficients of a polynomial expansion of the wave function.
//                The coefficients of the polynomials can be shown to decay exponentially for a well-resolved function.
//                Therefore, one can check the slope via a linear least-squares fit to (polynomial order, log(abs(polynomial coefficient)).
//                If the function is not well-resolved the polynomials decay more slowly and the slope is smaller.
//
// Note        :  1. This function is called once per patch group
//                2. The size of the array Var1D must be PS2 + 2
//                3. Assume a ghost size of 1
//
// Parameter   :  Var1D     : Array storing the input re & im
//                Cond      : Reference to floating point variable where density ratio is stored
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Prepare_for_Spectral_Criterion(const real *Var1D, real& Cond)
{
// set the stride to a small value to sample the wave function evenly
   const size_t Stride    = 1;
   const size_t GhostSize = 1;
   const size_t Size1D    = PS2 + 2 * GhostSize;
   const size_t MaxOrder  = 14;
   const size_t NField    = 3;
   const size_t NCoeff    = 2*NField;

   const real* Re1D = Var1D;
   const real* Im1D = Var1D + CUBE(Size1D);

   flag_spectral_float Order[MaxOrder];
   flag_spectral_float Coeff[NCoeff][MaxOrder]; // left and right coefficients for every field
   real Row[NField][Size1D];
   flag_spectral_float SlopeFront, SlopeBack, Intercept;

// initialise with large negative number
   Cond = -__FLT_MAX__;

// iterate over 3 dimensions and sample the physical 2D arrays with a stride
   for (size_t XYZ = 0; XYZ < 3; ++XYZ)
   for (size_t k=GhostSize; k<Size1D-GhostSize; k+=Stride)
   for (size_t j=GhostSize; j<Size1D-GhostSize; j+=Stride)
   {
//    read one column of data from 3D block
      for (size_t i = 0; i < Size1D; ++i) {
         size_t index;

         switch (XYZ)
         {
            case 0:
               index = IDX321(k, j, i, Size1D, Size1D);
               break;
            case 1:
               index = IDX321(k, i, j, Size1D, Size1D);
               break;
            case 2:
               index = IDX321(i, k, j, Size1D, Size1D);
               break;
         }

         Row[0][i] = Re1D[index];
         Row[1][i] = Im1D[index];
         Row[2][i] = SQR(Row[0][i]) + SQR(Row[1][i]);
      }

      for (int i = 0; i < MaxOrder; ++i)
      {
         for (int j = 0; j < NCoeff; ++j) {
            Coeff[j][i] = 0;
         }

//       Compute polynomial expansions of real and imaginary parts
         for (int t = 0; t < MaxOrder; t++) {
            for (int l = 0; l < NField; l++) {
               Coeff[l       ][i] += Flag_Spectral_Polynomials[i][t] * Row[l][t];                     // left boundary
               Coeff[l+NField][i] += Flag_Spectral_Polynomials[i][t] * Row[l][Size1D - MaxOrder + t]; // right boundary
            }
         } // t

//       Prepare linear fit to logarithm of polynomial coefficients
         Order[i] = log(i + 1);

         for (int j = 0; j < NCoeff; ++j) {
            Coeff[j][i] = log(abs(Coeff[j][i]) + 1e-16);
         }
      }

//    Find maximum slope to determine whether refinement is necessary
//    Large negative slopes indicate that wavefunction is well-resolved
      for (int j = 0; j < NCoeff; ++j) {
         const int n_lr = 8; // size of linear regression domain

//       Compute slope for first and last n_lr elements and take minimum float to
         Least_Squares_Regression(Order, Coeff[j], 0,                        MIN(n_lr , MaxOrder), &SlopeFront, &Intercept);
         Least_Squares_Regression(Order, Coeff[j], MAX(0, MaxOrder - n_lr ), MaxOrder,             &SlopeBack,  &Intercept);
         Cond = MAX(Cond, MIN(SlopeFront, SlopeBack));
      }
      printf("Cond: %f\n", Cond);

   } // XYZ, k,j
} // FUNCTION : Prepare_for_Spectral_Criterion


#endif // #if ( MODEL == ELBDM )
