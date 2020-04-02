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
    ui->btn_find_exons->setDisabled(true);
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
    if (this->sequence_path != ui->txt_UnknownDNAS_path->text()
            || this->dna_path != ui->txt_mRNA_path->text()
            || this->match != ui->spBox_match->value()
            || this->mismatch != ui->spBox_mismatch->value()
            || this->gap != ui->spBox_gap->value()
            || this->threshold != ui->txt_threshold->text())
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

        //1% progress
        QTextStream in_seq(&sequence_file);
        while (!in_seq.atEnd())
        {
            QString temp = in_seq.readLine();
            if (temp[0] == '>') continue;
            dna_sequence += temp;
        }
        sequence_file.close();

        QTextStream in_tpl(&dna_file);
        while (!in_tpl.atEnd())
        {
            QString temp = in_tpl.readLine();
            if (temp[0] == '>') continue;
            dna_template += temp;
        }

        //Progress diaglog
        QProgressDialog progress("Alignment in progress...", "Cancel", 0, dna_template.toStdString().size()*dna_sequence.toStdString().size()-1, this);
        progress.setMinimumDuration(100);
        progress.setWindowModality(Qt::WindowModal);

        clearTable();
        if (this->la) delete this->la;
        //98% progress
        this->la = new LocalAlignment(dna_template.toStdString(), dna_sequence.toStdString(),
                                      ui->spBox_match->value(), ui->spBox_mismatch->value(),
                                      ui->spBox_gap->value(), ui->txt_threshold->text().toInt(), progress);
        if (!progress.wasCanceled())
        {
            if (la->getExons().size() ==0)
            {
                QMessageBox::information(this, tr("Nothing found"), tr("There is no result!"));
                return;
            }
            else
            {
                this->sequence_path = ui->txt_UnknownDNAS_path->text();
                this->dna_path = ui->txt_mRNA_path->text();
                this->match = ui->spBox_match->value();
                this->mismatch = ui->spBox_mismatch->value();
                this->gap = ui->spBox_gap->value();
                this->threshold = ui->txt_threshold->text();

                qDebug() << "populate new data" << endl;
                populateTable(la->getExons());
                ui->btn_find_exons->setEnabled(true);
                this->reRunEC = true;
            }

        }
        else
        {
            qDebug() << "cancelled" << endl;
        }

    }
    else
    {
        if (this->tableValue != 1)
        {
            clearTable();
            populateTable(la->getExons());
            ui->btn_find_exons->setEnabled(true);
            qDebug() << "populate existing data" << endl;
        }
        else
        {
            qDebug() << "nothing todo" << endl;
        }
    }
    this->tableValue = 1;
}

void MainWindow::on_btn_find_exons_clicked()
{    
    if (!this->la || this->la->getExons().size() == 0)
    {
        QMessageBox::critical(this, tr("Error finding exons"), tr("The two sequences are not aligned yet"));
        return;
    }    

    if (this->reRunEC)
    {
        //Progress diaglog
        QProgressDialog progress("Exon chaining in progress...", "Cancel", 0, 100, this);
        progress.setMinimumDuration(100);
        progress.setWindowModality(Qt::WindowModal);

        clearTable();
        if (this->ec) delete this->ec;
        this->ec = new ExonChaining(this->la->getExons(), progress);
        progress.setValue(100);
        if (ec->get_intervals().size() == 0)
        {
            QMessageBox::information(this, tr("Nothing found"), tr("There is no result!"));
            return;
        }
        else
        {
            populateTable(ec->get_intervals());
            qDebug() << "populate new data" << endl;
        }
        this->reRunEC = false;
    }
    else
    {
        if (this->tableValue != 2)
        {
            clearTable();
            populateTable(ec->get_intervals());
            qDebug() << "populate existing data" << endl;
        }
        else
        {
            qDebug() << "nothing todo" << endl;
        }
    }
    this->tableValue = 2;
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
    unsigned int size = items.size();
    ui->lbl_count->setText("Total: " + QString::number(size));
    //Progress diaglog
    QProgressDialog progress("Data populating in progress...", "Cancel", 0, size-1, this);
    progress.setMinimumDuration(100);
    progress.setWindowModality(Qt::WindowModal);

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
        progress.setValue(i);
        ++i;
    }
}

void MainWindow::clearTable()
{
    ui->tbl_result->clearContents();
    ui->tbl_result->setRowCount(0);
    ui->lbl_count->setText("Total: 0");
    ui->ptxt_alignment->clear();
}

void MainWindow::on_txt_UnknownDNAS_path_textChanged(const QString &arg1)
{
    if (this->sequence_path != arg1)
    {
        ui->btn_find_exons->setDisabled(true);
    }
    else
    {
        ui->btn_find_exons->setEnabled(true);
    }
}

void MainWindow::on_txt_mRNA_path_textChanged(const QString &arg1)
{
    if (this->dna_path != arg1)
    {
        ui->btn_find_exons->setDisabled(true);
    }
    else
    {
        ui->btn_find_exons->setEnabled(true);
    }
}

void MainWindow::on_spBox_match_valueChanged(int arg1)
{
    if (this->match != arg1)
    {
        ui->btn_find_exons->setDisabled(true);
    }
    else
    {
        ui->btn_find_exons->setEnabled(true);
    }
}

void MainWindow::on_spBox_mismatch_valueChanged(int arg1)
{
    if (this->mismatch != arg1)
    {
        ui->btn_find_exons->setDisabled(true);
    }
    else
    {
        ui->btn_find_exons->setEnabled(true);
    }
}

void MainWindow::on_spBox_gap_valueChanged(int arg1)
{
    if (this->gap != arg1)
    {
        ui->btn_find_exons->setDisabled(true);
    }
    else
    {
        ui->btn_find_exons->setEnabled(true);
    }
}

void MainWindow::on_txt_threshold_textChanged(const QString &arg1)
{
    if (this->threshold != arg1)
    {
        ui->btn_find_exons->setDisabled(true);
    }
    else
    {
        ui->btn_find_exons->setEnabled(true);
    }
}
