# 원고 Word 파일을 PDF로 내보낸다(검토용). 사용: powershell -File export_pdf.ps1
$d = Split-Path -Parent $MyInvocation.MyCommand.Path
$w = New-Object -ComObject Word.Application; $w.Visible = $false
foreach ($n in @("review_copy_v2", "02_Supplementary_v2", "01_Main_Text_v2")) {
    $doc = $w.Documents.Open((Join-Path $d "$n.docx")); $doc.SaveAs([ref](Join-Path $d "$n.pdf"), [ref]17); $doc.Close()
}
$w.Quit()
