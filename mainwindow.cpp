#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->tbl_result->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->txt_threshold->setValidator(new QIntValidator(ui->txt_threshold));
    connect(ui->tbl_result, SIGNAL(cellDoubleClicked(int,int)), this, SLOT(tableItemDoubleClicked(int,int)));
}

MainWindow::~MainWindow()
{
    delete ui;
    delete la;
    delete ec;
}


void MainWindow::on_btn_browseUnknownDNA_file_clicked()
{
    QString filters = "All Files (*.*) ;; Text Files (*.txt) ;; FASTA Files (*.fasta)";
    QString file_name = QFileDialog::getOpenFileName(this, "Open Unknown DNA file", QDir::currentPath(), filters);
    ui->txt_UnknownDNAS_path->setText(file_name);
}

void MainWindow::on_btn_browseDNATemplate_file_clicked()
{
    QString filters = "All Files (*.*) ;; Text Files (*.txt) ;; FASTA Files (*.fasta)";
    QString file_name = QFileDialog::getOpenFileName(this, "Open mRNA Template file", QDir::currentPath(), filters);
    ui->txt_mRNA_path->setText(file_name);
}

void MainWindow::on_btn_align_clicked()
{
    if (this->reRunLA)
    {
        QString dna_template = "";
        QString dna_sequence = "";


        QFile sequence_file(ui->txt_UnknownDNAS_path->text());
        if (!sequence_file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            QMessageBox::critical(this, tr("Error opening file"), tr("Cannot open DNA sequence file"));
            return;
        }

        QFile dna_file(ui->txt_mRNA_path->text());
        if (!dna_file.open(QIODevice::ReadOnly | QIODevice::Text))
        {
            QMessageBox::critical(this, tr("Error opening file"), tr("Cannot open mRNA template file"));
            return;
        }

        QTextStream in_seq(&sequence_file);
        while (!in_seq.atEnd())
        {
            dna_sequence += in_seq.readLine();
        }
        sequence_file.close();

        QTextStream in_tpl(&dna_file);
        while (!in_tpl.atEnd())
        {
            dna_template += in_tpl.readLine();
        }
        if (this->la) delete this->la;
        this->la = new LocalAlignment(dna_template.toStdString(), dna_sequence.toStdString(),
                                      ui->spBox_match->text().toInt(), ui->spBox_mismatch->text().toInt(),
                                      ui->spBox_gap->text().toInt(), ui->txt_threshold->text().toInt());
        qDebug() << "rerun local alignment true";
        this->reRunEC = true;
    }
    populateTable(la->getExons());
    ui->ptxt_alignment->clear();
    this->reRunLA = false;

}

void MainWindow::on_btn_find_exons_clicked()
{
    if (!this->la)
    {
        QMessageBox::critical(this, tr("Error finding exons"), tr("The two sequences are not aligned yet"));
        return;
    }
    if (this->reRunEC)
    {
        if (this->ec) delete this->ec;
        this->ec = new ExonChaining(this->la->getExons());
        qDebug() << "rerun exon chaining true";
    }
    populateTable(ec->get_intervals());
    ui->ptxt_alignment->clear();
    this->reRunEC = false;
}

void MainWindow::tableItemDoubleClicked(int row, int col)
{
    if (col ==3)
    {

        unsigned int count = 0;
        QStringList start_s = ui->tbl_result->item(row, 0)->text().split('-');
        QStringList end_s = ui->tbl_result->item(row, 1)->text().split('-');
        Coordinate start(start_s[0].toInt(), start_s[1].toInt());
        Coordinate end(end_s[0].toInt(), end_s[1].toInt());
        std::list<std::array<std::string, 3>> alignments = this->la->print_alignment(start, end);
        /*QDialog dialog;
        dialog.setModal(true);
        dialog.exec();*/
        QString result = "";
        for (const std::array<std::string, 3> &alignment : alignments)
        {
            if (count > 0)
            {
                result += "\n\n";
                result +=  "------------------********new_alignment********------------------";
                result += "\n\n";
            }
            unsigned int char_breaks_count = 90;
            for (unsigned int i = 0; i < alignment[0].size(); i += char_breaks_count)
            {
                char_breaks_count = i+char_breaks_count > alignment[0].size() ? alignment[0].size() : char_breaks_count;
                if ( i >  0)
                {
                    result += '\n';
                    result +=  "------------------**break_line**------------------";
                    result += "\n\n";
                }
                result += QString::fromStdString(alignment[0].substr(i, char_breaks_count));
                result += '\n';
                result += QString::fromStdString(alignment[1].substr(i, char_breaks_count));
                result += '\n';
                result += QString::fromStdString(alignment[2].substr(i, char_breaks_count));
                result += '\n';
            }
            ++count;
        }
        ui->ptxt_alignment->setPlainText(result);
    }

}

void MainWindow::populateTable(const std::list<Interval_Coordinate> &items)
{
    ui->tbl_result->clearContents();
    unsigned int size = items.size();
    ui->tbl_result->setRowCount(size);
    unsigned int i = 0;
    for (const Interval_Coordinate &exon : items)
    {
        QString start = QString::number(exon.start.row) + "-" + QString::number(exon.start.col);
        QTableWidgetItem *start_item = new QTableWidgetItem(start);
        start_item->setTextAlignment(Qt::AlignCenter);
        start_item->setToolTip("template_position-sequence_position");
        ui->tbl_result->setItem(i,0,  start_item);

        QString end = QString::number(exon.end.row) + "-" + QString::number(exon.end.col);
        QTableWidgetItem *end_item = new QTableWidgetItem(end);
        end_item->setTextAlignment(Qt::AlignCenter);
        end_item->setToolTip("template_position-sequence_position");
        ui->tbl_result->setItem(i,1,  end_item);

        QTableWidgetItem *score_item = new QTableWidgetItem(QString::number(exon.score));
        score_item->setTextAlignment(Qt::AlignCenter);
        score_item->setToolTip("Alignment score");
        ui->tbl_result->setItem(i,2,  score_item);

        QTableWidgetItem *view_alignment = new QTableWidgetItem("View Alignment");
        QFont view_font;
        view_font.setUnderline(true);
        view_alignment->setFont(view_font);
        view_alignment->setTextAlignment(Qt::AlignCenter);
        view_alignment->setToolTip("Double click to view alignment");
        ui->tbl_result->setItem(i,3,  view_alignment);
        ++i;
    }
    ui->lbl_count->setText("Total: " + QString::number(size));
}

void MainWindow::on_txt_UnknownDNAS_path_textChanged(const QString &arg1)
{
    this->reRunLA = true;
    this->reRunEC = true;
}

void MainWindow::on_txt_mRNA_path_textChanged(const QString &arg1)
{
    this->reRunLA = true;
    this->reRunEC = true;
}

void MainWindow::on_spBox_match_valueChanged(int arg1)
{
    this->reRunLA = true;
    this->reRunEC = true;
}

void MainWindow::on_spBox_mismatch_valueChanged(int arg1)
{
    this->reRunLA = true;
    this->reRunEC = true;
}

void MainWindow::on_spBox_gap_valueChanged(int arg1)
{
    this->reRunLA = true;
    this->reRunEC = true;
}

void MainWindow::on_txt_threshold_textChanged(const QString &arg1)
{
    this->reRunLA = true;
    this->reRunEC = true;
}
