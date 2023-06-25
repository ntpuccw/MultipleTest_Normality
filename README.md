# MultipleTestNormality

Paper : <a href=http://jcsa.stat.org.tw/jcsa/data/vol61/V61N1-2.pdf>多變量常態的多重檢定實證研究</a><br>
 MATLAB codes: 以目錄分類

<ol>
 <li>src: 論文主程式</li>
 <ul>
  <li>prepaer_CV_from_tests.m: 計算論文式(10) 中的 C_n,h(u) 值，分別來自統計量 MK_lower, MK_upper, MS, BhS, BhL, Wmin(5)。計算結果儲存在目錄 new_data 以 Cnhu_p_ 為首的檔案。</li>
  <li>comp_cnu.m: 計算組合統計量在 0.0001:0.0001:0.3 每個位置的關鍵值  critical values。計算結果儲存在 new_data 目錄內以 Cnu_ 為首的檔案。</li>
  <li>compu_power: 根據不同的對立假設，計算三個多重檢定（MW, MBW, MMBB）在不同條件下的檢定力 power。 計算結果儲存在目錄 new_data 以 POWER_ 為首的檔案。</li>
 </ul>

 <li>presentations: 論文中的圖、表程式碼</li>
 <ul>
  <li>show_power.mlx:  繪製檢定力圖，如論文之圖 2 ~ 圖 6</li>
  <li>show_U_na_table.m: 計算論文表 2、表 3 的 U_n,a 值</li>
  <li>plot_Burr_Pareto_Logistic: 對立假說 Burr_Pareto_Logistic 分配的立體圖、等高線圖與樣本點的散佈</li>
  <li>plot_multi_Khintchine_GEP: 對立假說 multi_Khintchine_GEP 分配的立體圖、等高線圖與樣本點的散佈</li>
  <li>plot_multi_Khintchine_normal: 對立假說 Khintchine_normal 分配的立體圖、等高線圖與樣本點的散佈</li>
  <li>plot_multi_mmvn: 對立假說 mixed MVN 分配的立體圖、等高線圖與樣本點的散佈</li>
  <li>plot_multi_Pearson: 對立假說 multi_Pearson 分配的立體圖、等高線圖與樣本點的散佈</li>
  <li>plot_multi_SN: 對立假說 multi_Skewed-Normal 分配的立體圖、等高線圖與樣本點的散佈</li>
  <li>plot_multi_T: 對立假說 Multi_T 分配的立體圖、等高線圖與樣本點的散佈</li>
 </ul>
 <li>Tools: 論文使用的副程式：含自製與 open sources</li>
 <li>Miscell: 其他</li>
</ol>
