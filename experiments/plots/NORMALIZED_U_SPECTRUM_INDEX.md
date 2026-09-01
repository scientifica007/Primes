# Normalized-u experiments — SVG index

التقارير المرتبطة:

- [`../NORMALIZED_U_SPECTRUM_TRANSFER_STUDY.tex`](../NORMALIZED_U_SPECTRUM_TRANSFER_STUDY.tex)
- [`../NORMALIZED_U_POINT_MODEL_STUDY.tex`](../NORMALIZED_U_POINT_MODEL_STUDY.tex)

## 1. الطيف المرجعي في الإحداثي المحلي المطبع

![Normalized-u reference spectrum](normalized_u_reference_spectrum.svg)

- [`normalized_u_reference_spectrum.svg`](normalized_u_reference_spectrum.svg)

كل فترة بين مربعي أوليين متتاليين لها طول 1 في الإحداثي `u`. يبين الرسم أحد حلول NNLS للطيف المرجعي؛ القيم الفردية ليست فريدة بالضرورة، ولذلك لا تُفسر كقانون مستقل.

## 2. مقارنة النقل القديم بالنقل بعد التطبيع المحلي

![Normalized-u transfer improvement](normalized_u_transfer_improvement.svg)

- [`normalized_u_transfer_improvement.svg`](normalized_u_transfer_improvement.svg)

النقل في الإحداثي الفيزيائي باستخدام `L/H_median` كان يعطي متوسط kernel RMSE يقارب `2.60e-2`. بعد تحويل كل فترة إلى وحدة واحدة في `u` انخفض متوسط الخطأ على 17 لوحة غير مستخدمة في الملاءمة إلى نحو `4.88e-4` في النسخة uniform، أي تحسن بعامل يقارب 53.

## 3. نقل spectrum واحد من 10^7 إلى 10^12

![Normalized-u high-X kernel transfer](normalized_u_high_X_kernel_transfer.svg)

- [`normalized_u_high_X_kernel_transfer.svg`](normalized_u_high_X_kernel_transfer.svg)

الخط المتصل هو target kernel النظري في ثلاث لوحات قرب `2e12`، والمتقطع هو prediction الناتج من spectrum واحد ملائم قرب `2e7` مع تعديل amplitude اللوغاريتمي فقط ومن دون إعادة ملاءمة شكل الطيف.

## 4. إدخال spectrum المجمد إلى النموذج النقطي

![Normalized-u point model rho](normalized_u_point_model_full_rho.svg)

- [`normalized_u_point_model_full_rho.svg`](normalized_u_point_model_full_rho.svg)

يقارن منحنى الارتباط المرصود، short-interval theory، ومتوسط النموذج النقطي الذي يطبق fine spectrum المجمد مباشرة على ناجي random sieve. تمت معايرة معامل شدة واحد فقط من `R` على training set.

## 5. التشتت مقابل خطأ pair kernel

![Normalized-u point-model tradeoff](normalized_u_point_model_R_kernel_tradeoff.svg)

- [`normalized_u_point_model_R_kernel_tradeoff.svg`](normalized_u_point_model_R_kernel_tradeoff.svg)

يوضح أن normalized-u field يرفع التشتت إلى المجال المرصود ويحسن kernel الخام، لكنه لا يصل بعد إلى دقة interval pair-kernel bridge.

## 6. الفرق بين نجاح نقل النظرية ونجاح مولد النقاط

![Normalized-u theory vs point-model gap](normalized_u_theory_vs_point_model_gap.svg)

- [`normalized_u_theory_vs_point_model_gap.svg`](normalized_u_theory_vs_point_model_gap.svg)

التحويل إلى `u` خفض خطأ نقل fine spectrum النظري بنحو عامل 53، لكن داخل point generator بقي kernel RMSE قريبًا من physical local model. النتيجة تشير إلى أن التصحيح التالي يجب أن يستهدف residual covariance بعد طرح مساهمة random sieve، لا total covariance كاملة.
