#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "ExonChaining.h"
#include "LocalAlignment.h"
#include <QFileDialog>
#include <QDir>
#include <QMessageBox>
#include <QValidator>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_btn_browseUnknownDNA_file_clicked();
    void on_btn_browseDNATemplate_file_clicked();
    void on_btn_align_clicked();
    void tableItemDoubleClicked(int row, int col);
    void on_btn_find_exons_clicked();

    void on_txt_UnknownDNAS_path_textChanged(const QString &arg1);

    void on_txt_mRNA_path_textChanged(const QString &arg1);

    void on_spBox_match_valueChanged(int arg1);

    void on_spBox_mismatch_valueChanged(int arg1);

    void on_spBox_gap_valueChanged(int arg1);

    void on_txt_threshold_textChanged(const QString &arg1);

private:
    Ui::MainWindow *ui;
    LocalAlignment *la = nullptr;
    ExonChaining *ec = nullptr;
    bool reRunLA = true;
    bool reRunEC = true;
    void populateTable(const std::list<Interval_Coordinate> &items);
};
#endif // MAINWINDOW_H
